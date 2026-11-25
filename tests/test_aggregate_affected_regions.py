"""Verify affected-region cohort exports, source validation and safe writes."""

import contextlib
import csv
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "bin/3.1_T2Processing/aggregate_affected_regions.py"
spec = importlib.util.spec_from_file_location("aggregate_affected_regions", SCRIPT)
aggregate = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = aggregate
spec.loader.exec_module(aggregate)

HEADER = ["Label_id", "Brain region", "Affected region percentage [%]",
          "Region stroke volume [mm3]", "", "Stroke percent of brain [%]",
          "Total stroke volume [mm3]"]


class AffectedRegionExportTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)
        self.root = self.folder / "proc"
        self.output = self.folder / "exports" / "summary.csv"
        self.groups = self.folder / "groups.csv"
        self.groups.write_text("tDCS,Sham\nsub-NR1044,sub-B\n,sub-C\n", encoding="utf-8")
        self.translation = self.folder / "acronyms.txt"
        self.translation.write_text("22\t\tPTLp\n337 SSp-ll\n", encoding="utf-8")

    def sample(self, subject="sub-NR1044", session="ses-PT3", suffix="affectedRegions_Parental.csv",
               regions=(("L_Agranular_insular_area", "0.0818", "3.15"),), rows=None, header=None):
        directory = self.root / subject / session / "anat" / "affected_Regions"
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"{subject}_{session}_T2w_annoVolume+2000_rsfMRI_{suffix}"
        if rows is None:
            rows = [[index, name, percent, volume, "",
                     "2.73" if index == 1 else "", "14.4933" if index == 1 else ""]
                    for index, (name, volume, percent) in enumerate(regions, 1)]
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(HEADER if header is None else header)
            writer.writerows(rows)
        return path

    def export(self, **kwargs):
        options = dict(root=self.root, suffix="parental.csv", output=self.output)
        options.update(kwargs)
        summary = aggregate.export_cohort(**options)
        with self.output.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            return summary, reader.fieldnames, list(reader)

    def test_blocks_align_regions_and_sessions_and_preserve_quoted_names_and_precision(self):
        quoted = "L_Primary_somatosensory_area,_barrel_field"
        first = self.sample(regions=((quoted, "0.8925", "30.83"),
                                     ("L_Agranular_insular_area", "0.0818", "3.15")))
        original = first.read_bytes()
        self.sample(session="ses-PT0", regions=(("R_Region", "0.123456789123456789", "0.0001"),
                                               (quoted, "0.1111", "1.01")))
        self.sample("sub-B", regions=((quoted, "0.2000", "2.00"),))
        summary, header, rows = self.export(group_csv=self.groups)
        names = ["L_Agranular_insular_area", quoted, "R_Region"]
        self.assertEqual(header, ["subject", "session", "group", *aggregate.TOTALS]
                         + [f"{name} [mm3]" for name in names]
                         + [f"{name} [%]" for name in names])
        self.assertEqual((summary.sessions, summary.subjects, summary.regions), (3, 2, 3))
        self.assertEqual([(row["subject"], row["session"]) for row in rows],
                         [("sub-B", "ses-PT3"), ("sub-NR1044", "ses-PT0"), ("sub-NR1044", "ses-PT3")])
        self.assertEqual([row["group"] for row in rows], ["Sham", "tDCS", "tDCS"])
        self.assertEqual(rows[1]["R_Region [mm3]"], "0.123456789123456789")
        self.assertEqual(rows[1]["L_Agranular_insular_area [mm3]"], "0")
        self.assertEqual(rows[1]["L_Agranular_insular_area [%]"], "0")
        self.assertEqual(rows[2][f"{quoted} [mm3]"], "0.8925")
        self.assertEqual(rows[2][f"{quoted} [%]"], "30.83")
        self.assertEqual(rows[2]["Total stroke volume [mm3]"], "14.4933")
        self.assertEqual(rows[2]["Stroke percent of brain [%]"], "2.73")
        self.assertEqual(first.read_bytes(), original)

    def test_no_group_argument_omits_group_column(self):
        self.sample()
        _, header, rows = self.export()
        self.assertEqual(header[:4], ["subject", "session", *aggregate.TOTALS])
        self.assertNotIn("group", header)
        self.assertNotIn("group", rows[0])

    def test_unassigned_subject_is_retained_with_warning_once_for_all_sessions(self):
        self.sample("sub-Unknown", "ses-PT0")
        self.sample("sub-Unknown", "ses-PT3")
        warnings = io.StringIO()
        with contextlib.redirect_stderr(warnings):
            _, _, rows = self.export(group_csv=self.groups)
        self.assertEqual([row["group"] for row in rows], ["", ""])
        self.assertEqual(warnings.getvalue().count("No group assignment for sub-Unknown"), 1)

    def test_duplicate_group_subject_is_rejected_even_within_one_group(self):
        self.sample()
        for text in ("tDCS,Sham\nsub-NR1044,sub-NR1044\n", "tDCS\nsub-NR1044\nsub-NR1044\n"):
            with self.subTest(text=text):
                self.groups.write_text(text, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "occurs more than once"):
                    self.export(group_csv=self.groups)
                self.assertFalse(self.output.exists())

    def test_group_csv_handles_bom_whitespace_and_uneven_columns(self):
        self.sample()
        self.sample("sub-C")
        self.groups.write_text("\ufefftDCS;Sham\n sub-NR1044 ; sub-B\n;sub-C\n", encoding="utf-8")
        _, _, rows = self.export(group_csv=self.groups)
        self.assertEqual([row["group"] for row in rows], ["Sham", "tDCS"])

    def test_suffix_is_literal_case_insensitive_and_selects_only_matching_atlas(self):
        self.sample(suffix="affectedRegions_Parental.CSV")
        self.sample(suffix="affectedRegions.csv", regions=(("Detailed", "1.00", "5.00"),))
        directory = self.root / "sub-NR1044/ses-PT3/anat/affected_Regions"
        (directory / "unrelated_parental.csv").write_text("not an affectedRegions CSV", encoding="utf-8")
        _, header, _ = self.export(suffix="PARENTAL.csv")
        self.assertIn("L_Agranular_insular_area [mm3]", header)
        self.assertNotIn("Detailed [mm3]", header)
        _, header, _ = self.export(suffix="affectedRegions.csv")
        self.assertIn("Detailed [mm3]", header)
        self.assertNotIn("L_Agranular_insular_area [mm3]", header)

    def test_missing_sessions_and_anat_folders_are_skipped_with_warning(self):
        self.sample()
        (self.root / "sub-B/ses-PT0").mkdir(parents=True)
        self.sample("sub-C", suffix="affectedRegions.csv")
        warnings = io.StringIO()
        with contextlib.redirect_stderr(warnings):
            summary, _, _ = self.export()
        self.assertEqual(summary.sessions, 1)
        self.assertIn("Skipping sub-B/ses-PT0", warnings.getvalue())
        self.assertIn("Skipping sub-C/ses-PT3", warnings.getvalue())

    def test_multiple_matches_are_rejected_with_filenames(self):
        path = self.sample()
        extra = path.with_name("other_affectedRegions_Parental.csv")
        extra.write_bytes(path.read_bytes())
        with self.assertRaisesRegex(ValueError, "Multiple matching CSVs") as caught:
            self.export()
        self.assertIn(extra.name, str(caught.exception))
        self.assertFalse(self.output.exists())

    def test_no_matches_and_invalid_root_are_errors(self):
        with self.assertRaisesRegex(ValueError, "Input folder does not exist"):
            self.export()
        self.root.mkdir()
        with self.assertRaisesRegex(ValueError, "No matching affectedRegions CSVs"):
            self.export()

    def test_invalid_sources_preserve_existing_output_and_report_path(self):
        good = [1, "L_A", "3.15", "0.0818", "", "2.73", "14.4933"]
        cases = [
            (HEADER[:-1], [good[:-1]], "Missing required columns"),
            (HEADER, [], "No data rows"),
            (HEADER, [good, good], "duplicate brain region"),
            (HEADER, [[1, "", "3.15", "0.0818", "", "2.73", "14.4933"]], "empty brain region"),
            (HEADER, [[1, "L_A", "", "0.0818", "", "2.73", "14.4933"]], "expected a number"),
            (HEADER, [[1, "L_A", "NaN", "0.0818", "", "2.73", "14.4933"]], "finite number"),
            (HEADER, [[1, "L_A", "3.15", "inf", "", "2.73", "14.4933"]], "finite number"),
            (HEADER, [[1, "L_A", "3.15", "0.0818", "", "", ""]], "Missing total stroke values"),
            (HEADER, [good, [2, "L_B", "1.00", "0.1000", "", "2.74", "14.4933"]], "conflicting values"),
            (HEADER, [good, [2, "L_B", "1.00", "0.1000", "", "2.73", "14.5000"]], "conflicting values"),
            (HEADER, [[1, "L_A", "3.15", "0.0818", "", "NaN", "14.4933"]], "finite number"),
            (HEADER, [good[:-1]], "column count"),
            (HEADER, [good + ["extra"]], "column count"),
            (HEADER + ["Brain region"], [good + ["L_A"]], "Duplicate column names"),
        ]
        self.output.parent.mkdir(parents=True)
        self.output.write_bytes(b"existing output\n")
        for header, rows, message in cases:
            with self.subTest(message=message, rows=rows):
                path = self.sample(rows=rows, header=header)
                original = path.read_bytes()
                with self.assertRaisesRegex(ValueError, message) as caught:
                    self.export()
                self.assertIn(str(path), str(caught.exception))
                self.assertEqual(self.output.read_bytes(), b"existing output\n")
                self.assertEqual(path.read_bytes(), original)

    def test_totals_can_be_on_later_rows_or_repeated_without_rounding(self):
        self.sample(rows=[[1, "L_A", "1.00", "0.1111", "", "", ""],
                          [2, "L_B", "2.00", "0.2222", "", "2.730", "14.493300"],
                          [3, "L_C", "3.00", "0.3333", "", "2.73", "14.4933"]])
        _, _, rows = self.export()
        self.assertEqual(rows[0]["Total stroke volume [mm3]"], "14.493300")
        self.assertEqual(rows[0]["Stroke percent of brain [%]"], "2.730")

    def test_malformed_quoted_csv_reports_source_path(self):
        path = self.sample()
        with path.open("a", encoding="utf-8") as stream:
            stream.write('2,"unclosed region')
        with self.assertRaisesRegex(ValueError, "unexpected end of data") as caught:
            self.export()
        self.assertIn(str(path), str(caught.exception))

    def test_output_cannot_overwrite_source_or_group_csv(self):
        path = self.sample()
        for output in (path, self.groups):
            with self.subTest(output=output):
                original = output.read_bytes()
                with self.assertRaisesRegex(ValueError, "must not overwrite an input"):
                    self.export(output=output, group_csv=self.groups)
                self.assertEqual(output.read_bytes(), original)

    def test_failed_write_preserves_existing_output_and_removes_temporary_file(self):
        self.sample()
        self.output.parent.mkdir(parents=True)
        self.output.write_bytes(b"existing output\n")
        for operation in ("replace", "writer"):
            with self.subTest(operation=operation):
                target = aggregate.os if operation == "replace" else aggregate.csv
                with mock.patch.object(target, operation, side_effect=OSError("write failed")):
                    with self.assertRaisesRegex(OSError, "write failed"):
                        self.export()
                self.assertEqual(self.output.read_bytes(), b"existing output\n")
                self.assertEqual(list(self.output.parent.glob(".affected_regions_*.csv")), [])

    def test_cli_exports_and_reports_summary(self):
        self.sample()
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            aggregate.main(["-i", str(self.root), "-a", "parental.csv", "-o", str(self.output),
                            "--group", str(self.groups)])
        self.assertTrue(self.output.is_file())
        self.assertIn("Saved 1 sessions from 1 subjects, 1 regions", stdout.getvalue())

    def test_cli_invalid_input_has_readable_error_without_traceback(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), self.assertRaises(SystemExit) as caught:
            aggregate.main(["-i", str(self.root), "-a", "parental.csv", "-o", str(self.output)])
        self.assertEqual(caught.exception.code, 2)
        self.assertIn("Input folder does not exist", stderr.getvalue())
        self.assertNotIn("Traceback", stderr.getvalue())

    def test_translation_expands_both_hemispheres_and_preserves_case_and_source(self):
        self.translation.write_text("\ufeff22\t\tPTLp\n\n 337   SSp-ll  \n", encoding="utf-8")
        original = self.translation.read_bytes()
        self.assertEqual(aggregate.read_translation(self.translation),
                         {22: "L_PTLp", 2022: "R_PTLp", 337: "L_SSp-ll", 2337: "R_SSp-ll"})
        self.assertEqual(self.translation.read_bytes(), original)

    def test_translated_blocks_use_ids_align_samples_and_preserve_values(self):
        source = self.sample(rows=[
            [337, "AAA unrelated name", "23.54", "0.2266", "", "2.73", "14.4933"],
            [2022, "Identical original name", "5.00", "0.123456789123456789", "", "", ""],
            [22, "Identical original name", "6.00", "0.5000", "", "", ""],
        ])
        self.sample(session="ses-PT0", rows=[
            [337, "Different name for the same ID", "0.01", "0.0001", "", "2.73", "14.4933"],
            [2337, "Contradictory L_ name", "2.50", "0.4321", "", "", ""],
        ])
        original_source = source.read_bytes()
        original_translation = self.translation.read_bytes()
        summary, header, rows = self.export(translation_txt=self.translation, group_csv=self.groups)
        names = ["L_PTLp", "L_SSp-ll", "R_PTLp", "R_SSp-ll"]
        self.assertEqual(header, ["subject", "session", "group", *aggregate.TOTALS]
                         + [f"{name} [mm3]" for name in names]
                         + [f"{name} [%]" for name in names])
        self.assertEqual((summary.sessions, summary.regions), (2, 4))
        self.assertEqual(rows[0]["L_PTLp [mm3]"], "0")
        self.assertEqual(rows[0]["L_PTLp [%]"], "0")
        self.assertEqual(rows[0]["L_SSp-ll [mm3]"], "0.0001")
        self.assertEqual(rows[0]["R_SSp-ll [%]"], "2.50")
        self.assertEqual(rows[1]["L_PTLp [mm3]"], "0.5000")
        self.assertEqual(rows[1]["L_SSp-ll [%]"], "23.54")
        self.assertEqual(rows[1]["R_PTLp [mm3]"], "0.123456789123456789")
        self.assertEqual(rows[1]["R_SSp-ll [mm3]"], "0")
        self.assertEqual([row["group"] for row in rows], ["tDCS", "tDCS"])
        self.assertEqual(rows[1]["Total stroke volume [mm3]"], "14.4933")
        self.assertEqual(rows[1]["Stroke percent of brain [%]"], "2.73")
        self.assertEqual(source.read_bytes(), original_source)
        self.assertEqual(self.translation.read_bytes(), original_translation)

    def test_bad_translation_tables_report_location_and_preserve_output(self):
        self.sample()
        self.output.parent.mkdir(parents=True)
        self.output.write_bytes(b"existing output\n")
        cases = [
            ("\n", "contains no label ID/acronym pairs", None),
            ("22\n", "expected a label ID and acronym", 1),
            ("bad PTLp\n", "invalid integer label ID", 1),
            ("22 PTLp\n22 PTLp\n", "duplicate label ID 22", 2),
            ("22 PTLp\n22 ACA\n", "duplicate label ID 22", 2),
            ("22 PTLp\n31 PTLp\n", "ambiguous acronym", 2),
            ("22 PTLp\n2022 ACA\n", "ID collision at 2022", 2),
            ("2022 ACA\n22 PTLp\n", "ID collision at 2022", 2),
        ]
        for text, message, line in cases:
            with self.subTest(text=text):
                self.translation.write_text(text, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message) as caught:
                    self.export(translation_txt=self.translation)
                self.assertIn(str(self.translation), str(caught.exception))
                if line is not None:
                    self.assertIn(f"Line {line}", str(caught.exception))
                self.assertEqual(self.output.read_bytes(), b"existing output\n")
                self.assertEqual(self.translation.read_text(encoding="utf-8"), text)

    def test_invalid_or_missing_translation_ids_report_source_and_preserve_output(self):
        self.output.parent.mkdir(parents=True)
        self.output.write_bytes(b"existing output\n")
        for label_id, message in (("", "invalid integer Label_id"), ("337.0", "invalid integer Label_id"),
                                  ("bad", "invalid integer Label_id"), (999, "no translation for Label_id 999")):
            with self.subTest(label_id=label_id):
                path = self.sample(rows=[[label_id, "L_Region", "3.15", "0.0818", "", "2.73", "14.4933"]])
                original = path.read_bytes()
                with self.assertRaisesRegex(ValueError, message) as caught:
                    self.export(translation_txt=self.translation)
                self.assertIn(str(path), str(caught.exception))
                self.assertIn("Line 2", str(caught.exception))
                self.assertEqual(path.read_bytes(), original)
                self.assertEqual(self.output.read_bytes(), b"existing output\n")

    def test_label_id_column_is_required_only_with_translation(self):
        path = self.sample(header=HEADER[1:], rows=[["L_A", "3.15", "0.0818", "", "2.73", "14.4933"]])
        _, header, _ = self.export()
        self.assertIn("L_A [mm3]", header)
        original_output = self.output.read_bytes()
        with self.assertRaisesRegex(ValueError, "Missing required columns: Label_id") as caught:
            self.export(translation_txt=self.translation)
        self.assertIn(str(path), str(caught.exception))
        self.assertEqual(self.output.read_bytes(), original_output)

    def test_duplicate_source_ids_cannot_overwrite_translated_values(self):
        self.sample(rows=[[22, "First name", "1.00", "0.1000", "", "2.73", "14.4933"],
                          [22, "Different name", "2.00", "0.2000", "", "", ""]])
        with self.assertRaisesRegex(ValueError, "duplicate Label_id 22"):
            self.export(translation_txt=self.translation)
        self.assertFalse(self.output.exists())

    def test_output_cannot_overwrite_translation_txt(self):
        self.sample(rows=[[22, "L_Region", "3.15", "0.0818", "", "2.73", "14.4933"]])
        original = self.translation.read_bytes()
        with self.assertRaisesRegex(ValueError, "must not overwrite an input file"):
            self.export(translation_txt=self.translation, output=self.translation)
        self.assertEqual(self.translation.read_bytes(), original)

    def test_cli_accepts_both_translation_option_names(self):
        self.sample(rows=[[2022, "R_Region", "3.15", "0.0818", "", "2.73", "14.4933"]])
        for option in ("-t", "--translationTXT"):
            with self.subTest(option=option), contextlib.redirect_stdout(io.StringIO()):
                aggregate.main(["-i", str(self.root), "-a", "parental.csv", "-o", str(self.output),
                                option, str(self.translation)])
                with self.output.open(encoding="utf-8", newline="") as stream:
                    reader = csv.DictReader(stream)
                    rows = list(reader)
                self.assertIn("R_PTLp [mm3]", reader.fieldnames)
                self.assertEqual(rows[0]["R_PTLp [%]"], "3.15")


if __name__ == "__main__":
    unittest.main()
