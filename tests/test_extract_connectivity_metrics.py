"""Verify cohort selection, group mapping and regional CSV alignment."""

import contextlib
import csv
import importlib.util
import io
from pathlib import Path
import stat
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
from scipy.io import savemat


SCRIPT = Path(__file__).resolve().parents[1] / "bin/3.2_DTIConnectivity/extract_connectivity_metrics.py"
spec = importlib.util.spec_from_file_location("connectivity_metrics", SCRIPT)
extract = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = extract
spec.loader.exec_module(extract)


class ConnectivityMetricsTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)
        self.root = self.folder / "proc"
        self.groups = self.folder / "groups.csv"
        self.groups.write_text("tDCS,Sham\nsub-A,sub-B\n", encoding="utf-8")
        self.output = self.folder / "result.csv"

    def sample(self, subject, session="ses-PT0", names=("L_Primary_motor_area", "R_Primary_motor_area"),
               counts=(7, 5), matrix=None, atlas="parental", connectivity="pass", converted=False,
               resampled=False, metrics=None):
        folder = self.root / subject / session / "dwi" / "connectivity"
        folder.mkdir(parents=True, exist_ok=True)
        marker = "AnnoSplit_parental" if atlas == "parental" else "AnnoSplit"
        if resampled:
            marker += "_resampled"
        path = folder / f"{subject}_{session}_{marker}.all.{connectivity}.connectivity.mat"
        if matrix is None:
            matrix = np.array([[0, 3], [3, 0]], dtype=np.float32)
        data = {"number of tracts t2r": np.array(counts, dtype=np.float32)[:, None],
                "number of tracts r2r": np.asarray(matrix),
                "name": np.frombuffer(("\n".join(names) + "\n\0").encode(), dtype=np.uint8)[None, :]}
        if metrics:
            data.update(metrics)
        if converted:
            data = {extract.matlab_identifier(key): value for key, value in data.items()}
            data["region_names"] = np.asarray(names, dtype=object)[:, None]
        savemat(path, data, format="5" if converted else "4", do_compression=converted)
        return path

    def export(self, **kwargs):
        options = dict(atlas="parental", connectivity="pass", level="t2r")
        options.update(kwargs)
        summary = extract.export_cohort(self.root, self.groups, self.output, **options)
        with self.output.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            return summary, reader.fieldnames, list(reader)

    def test_t2r_aligns_names_groups_and_sessions_across_mat_formats(self):
        first = self.sample("sub-A")
        original = first.read_bytes()
        self.sample("sub-A", "ses-PT3", names=("R_Primary_motor_area", "L_Primary_motor_area"),
                    counts=(2, 9), converted=True)
        self.sample("sub-B", names=("L_Primary_motor_area", "L_corticospinal_tract"), counts=(0, 6))
        # Wrong atlas/type files must not be selected.
        self.sample("sub-A", atlas="detailed", counts=(99, 99))
        self.sample("sub-A", connectivity="end", counts=(99, 99))
        summary, header, rows = self.export()
        self.assertEqual((summary.sessions, summary.subjects, summary.columns), (3, 2, 6))
        self.assertEqual(header[:3], ["subject", "session", "group"])
        self.assertEqual(rows[0]["group"], "tDCS")
        self.assertEqual(rows[1]["L_Primary_motor_area"], "9")
        self.assertEqual(rows[1]["R_Primary_motor_area"], "2")
        self.assertEqual(rows[2]["group"], "Sham")
        self.assertEqual(rows[2]["R_Primary_motor_area"], "")
        self.assertEqual(rows[2]["L_Primary_motor_area"], "0")
        self.assertEqual(first.read_bytes(), original)

    def test_r2r_outputs_unique_pairs_without_diagonal_and_keeps_missing_empty(self):
        self.sample("sub-A")
        self.sample("sub-B", names=("R_Primary_motor_area", "L_corticospinal_tract"),
                    matrix=[[0, 4], [4, 0]], converted=True)
        summary, header, rows = self.export(level="r2r")
        self.assertEqual(summary.columns, 6)
        self.assertEqual(len(header[3:]), len(set(header[3:])))
        self.assertEqual(rows[0]["L_Primary_motor_area_x_R_Primary_motor_area"], "3")
        self.assertEqual(rows[1]["L_corticospinal_tract_x_R_Primary_motor_area"], "4")
        self.assertEqual(rows[1]["L_Primary_motor_area_x_R_Primary_motor_area"], "")
        for column in header[3:]:
            left, right = column.split("_x_")
            self.assertNotEqual(left, right)
            self.assertNotIn(right + "_x_" + left, header)

    def test_pair_accepts_exact_motor_and_corticospinal_names(self):
        self.sample("sub-A", names=("L_Primary_motor_area", "L_corticospinal_tract"))
        summary, header, rows = self.export(level="r2r", region1="L_Primary_motor_area",
                                            region2="L_corticospinal_tract")
        self.assertEqual(summary.columns, 4)
        self.assertEqual(header[3], "L_Primary_motor_area_x_L_corticospinal_tract")
        self.assertEqual(rows[0][header[3]], "3")
        _, reverse_header, reverse_rows = self.export(level="r2r", region1="L_corticospinal_tract",
                                                       region2="L_Primary_motor_area")
        self.assertEqual(reverse_header, header)
        self.assertEqual(reverse_rows, rows)

    def test_acronyms_and_inexact_names_are_rejected(self):
        self.sample("sub-A", names=("L_Primary_motor_area", "L_corticospinal_tract"))
        for query in ("L_MOp", "L_MoP", "L_cst", "l_primary_motor_area"):
            with self.subTest(query=query), self.assertRaisesRegex(ValueError, "full region name"):
                self.export(level="r2r", region1=query, region2="L_corticospinal_tract")
            self.assertFalse(self.output.exists())

    def test_detailed_end_selection_and_arbitrary_metric_preserve_precision(self):
        value = .12345678912345678
        names = ("L_Primary_motor_area,_Layer_5", "R_Primary_motor_area,_Layer_5")
        self.sample("sub-A", names=names, atlas="detailed", connectivity="end", resampled=True,
                    metrics={"mean length(mm) r2r": np.array([[0, value], [value, 0]])})
        self.sample("sub-B", names=names, atlas="detailed", connectivity="end", converted=True,
                    metrics={"mean length(mm) r2r": np.array([[0, value], [value, 0]])})
        _, header, rows = self.export(atlas="detailed", connectivity="end", level="r2r", metric="mean_length_mm",
                                     region1=names[0], region2=names[1])
        for row in rows:
            self.assertEqual(float(row[header[3]]), value)

    def test_csv_quotes_region_names_containing_commas(self):
        names = ('L_motor, "layer"', "R_motor ü")
        self.sample("sub-A", names=names)
        _, _, rows = self.export(level="r2r", region1=names[0], region2=names[1])
        self.assertEqual(rows[0]["_x_".join(names)], "3")

    def test_semicolon_bom_and_unequal_group_lengths(self):
        self.groups.write_text("\ufefftDCS;Sham\n sub-A ;sub-B\n;sub-C\n", encoding="utf-8")
        self.assertEqual(extract.read_groups(self.groups),
                         {"sub-A": "tDCS", "sub-B": "Sham", "sub-C": "Sham"})

    def test_missing_and_duplicate_groups_fail_before_export(self):
        self.sample("sub-A")
        for text, message in [("tDCS,Sham\nsub-A,sub-A\n", "more than once"),
                              ("tDCS,Sham\n,sub-B\n", "missing from group")]:
            with self.subTest(text=text):
                self.groups.write_text(text)
                with self.assertRaisesRegex(ValueError, message):
                    self.export()
                self.assertFalse(self.output.exists())

    def test_invalid_pair_arguments_and_regions_fail_before_export(self):
        self.sample("sub-A")
        cases = [({"level": "r2r", "region1": "L_Primary_motor_area"}, "both"),
                 ({"region1": "L_Primary_motor_area", "region2": "R_Primary_motor_area"}, "only"),
                 ({"level": "r2r", "region1": "L_Primary_motor_area", "region2": "L_Primary_motor_area"}, "different"),
                 ({"level": "r2r", "region1": "unknown", "region2": "R_Primary_motor_area"}, "absent")]
        for options, message in cases:
            with self.subTest(options=options):
                with self.assertRaisesRegex(ValueError, message):
                    self.export(**options)
                self.assertFalse(self.output.exists())

    def test_missing_mat_is_warned_and_skipped_for_both_levels(self):
        path = self.sample("sub-A")
        original = path.read_bytes()
        missing = self.root / "sub-B" / "ses-PT3" / "dwi"
        missing.mkdir(parents=True)
        for level in ("t2r", "r2r"):
            with self.subTest(level=level), contextlib.redirect_stderr(io.StringIO()) as stderr:
                summary, _, rows = self.export(level=level)
            self.assertEqual((summary.sessions, summary.subjects), (1, 1))
            self.assertEqual([row["subject"] for row in rows], ["sub-A"])
            self.assertIn("WARNING: Skipping sub-B/ses-PT3", stderr.getvalue())
            self.assertIn("parental/pass", stderr.getvalue())
            self.assertIn(str(missing / "connectivity"), stderr.getvalue())
        self.assertEqual(path.read_bytes(), original)

    def test_all_missing_mat_files_do_not_create_or_replace_output(self):
        # A MAT file for another connectivity type must count as missing.
        self.sample("sub-A", connectivity="end")
        for existing in (False, True):
            with self.subTest(existing=existing):
                if existing:
                    self.output.write_text("existing result\n")
                with contextlib.redirect_stderr(io.StringIO()) as stderr:
                    with self.assertRaisesRegex(ValueError, "No sessions with a matching parental/pass"):
                        self.export()
                self.assertIn("WARNING: Skipping sub-A/ses-PT0", stderr.getvalue())
                if existing:
                    self.assertEqual(self.output.read_text(), "existing result\n")
                else:
                    self.assertFalse(self.output.exists())
                self.assertEqual(list(self.folder.glob(".connectivity_*.csv")), [])

    def test_duplicate_mat_and_unknown_metric_are_reported(self):
        path = self.sample("sub-A")
        with self.assertRaisesRegex(ValueError, "Available metrics"):
            self.export(metric="not_a_metric")
        duplicate = path.with_name("duplicate_" + path.name)
        duplicate.write_bytes(path.read_bytes())
        with self.assertRaisesRegex(ValueError, "found 2"):
            self.export()
        self.assertFalse(self.output.exists())

    def test_corrupt_mat_still_aborts_without_replacing_output(self):
        path = self.sample("sub-A")
        corrupt = b"not a MATLAB file" * 16
        path.write_bytes(corrupt)
        self.output.write_text("existing result\n")
        with self.assertRaises(ValueError):
            self.export()
        self.assertEqual(self.output.read_text(), "existing result\n")
        self.assertEqual(path.read_bytes(), corrupt)

    def test_wrong_shape_and_duplicate_region_names_fail(self):
        path = self.sample("sub-A", matrix=np.zeros((3, 3)))
        with self.assertRaisesRegex(ValueError, "shape"):
            self.export(level="r2r")
        path.unlink()
        self.sample("sub-A", names=("L_Primary_motor_area", "L_Primary_motor_area"))
        with self.assertRaisesRegex(ValueError, "duplicate region"):
            self.export()

    def test_asymmetric_matrix_does_not_replace_output_and_temp_is_removed(self):
        self.sample("sub-A", matrix=[[0, 3], [4, 0]])
        self.output.write_text("existing result\n")
        with self.assertRaisesRegex(ValueError, "asymmetric"):
            self.export(level="r2r")
        self.assertEqual(self.output.read_text(), "existing result\n")
        self.assertEqual(list(self.folder.glob(".connectivity_*.csv")), [])

    def test_group_csv_cannot_be_overwritten(self):
        self.sample("sub-A")
        original = self.groups.read_bytes()
        with self.assertRaisesRegex(ValueError, "overwrite"):
            extract.export_cohort(self.root, self.groups, self.groups, "parental", "pass", "t2r")
        self.assertEqual(self.groups.read_bytes(), original)

    def test_new_and_replaced_csv_use_normal_permissions_under_umask(self):
        path = self.sample("sub-A")
        original = path.read_bytes()
        for mask, expected_mode in ((0o022, 0o644), (0o002, 0o664), (0o077, 0o600)):
            for existing in (False, True):
                with self.subTest(mask=oct(mask), existing=existing):
                    self.output.unlink(missing_ok=True)
                    if existing:
                        self.output.write_text("existing result\n")
                        self.output.chmod(0o600)
                    previous_umask = extract.os.umask(mask)
                    try:
                        _, _, rows = self.export()
                    finally:
                        extract.os.umask(previous_umask)
                    self.assertEqual(stat.S_IMODE(self.output.stat().st_mode), expected_mode)
                    self.assertEqual(rows[0]["L_Primary_motor_area"], "7")
                    self.assertEqual(path.read_bytes(), original)

    def test_replace_failure_preserves_output_and_removes_temporary_file(self):
        self.sample("sub-A")
        self.output.write_text("existing result\n")
        self.output.chmod(0o600)
        with mock.patch.object(extract.os, "replace", side_effect=PermissionError("replace denied")):
            with self.assertRaisesRegex(PermissionError, "replace denied"):
                self.export()
        self.assertEqual(self.output.read_text(), "existing result\n")
        self.assertEqual(stat.S_IMODE(self.output.stat().st_mode), 0o600)
        self.assertEqual(list(self.folder.glob(".connectivity_*.csv")), [])

    def test_write_failure_preserves_output_and_removes_temporary_file(self):
        self.sample("sub-A")
        self.output.write_text("existing result\n")
        with mock.patch.object(extract.csv, "writer") as writer:
            writer.return_value.writerow.side_effect = OSError("write failed")
            with self.assertRaisesRegex(OSError, "write failed"):
                self.export()
        self.assertEqual(self.output.read_text(), "existing result\n")
        self.assertEqual(list(self.folder.glob(".connectivity_*.csv")), [])

    def test_temporary_name_collision_preserves_existing_file_and_retries(self):
        self.sample("sub-A")
        first, second = extract.uuid.UUID(int=1), extract.uuid.UUID(int=2)
        collision = self.folder / f".connectivity_{first.hex}.csv"
        collision.write_text("another export\n")
        with mock.patch.object(extract.uuid, "uuid4", side_effect=[first, second]):
            summary, _, rows = self.export()
        self.assertEqual(summary.sessions, 1)
        self.assertEqual(rows[0]["L_Primary_motor_area"], "7")
        self.assertEqual(collision.read_text(), "another export\n")
        self.assertEqual(list(self.folder.glob(".connectivity_*.csv")), [collision])

    def test_cli_writes_the_requested_pair(self):
        self.sample("sub-A")
        (self.root / "sub-B" / "ses-PT0" / "dwi").mkdir(parents=True)
        with contextlib.redirect_stdout(io.StringIO()) as output, \
                contextlib.redirect_stderr(io.StringIO()) as stderr:
            extract.main(["--input", str(self.root), "--group", str(self.groups),
                          "--output", str(self.output), "--atlas", "parental", "--connectivity", "pass",
                          "--level", "r2r", "--region1", "L_Primary_motor_area",
                          "--region2", "R_Primary_motor_area"])
        self.assertIn("1 sessions from 1 subjects, 4 columns", output.getvalue())
        self.assertIn("WARNING: Skipping sub-B/ses-PT0", stderr.getvalue())
        with self.output.open(newline="") as stream:
            rows = list(csv.reader(stream))
        self.assertEqual(rows[1], ["sub-A", "ses-PT0", "tDCS", "3"])

    def test_cli_without_group_csv_omits_group_column_for_both_levels(self):
        self.sample("sub-A")
        self.groups.unlink()
        for level in ("t2r", "r2r"):
            with self.subTest(level=level), contextlib.redirect_stdout(io.StringIO()) as stdout:
                extract.main(["--input", str(self.root), "--output", str(self.output),
                              "--atlas", "parental", "--connectivity", "pass", "--level", level])
                with self.output.open(newline="") as stream:
                    reader = csv.DictReader(stream)
                    self.assertEqual(reader.fieldnames[:2], ["subject", "session"])
                    self.assertNotIn("group", reader.fieldnames)
                    self.assertEqual(len(reader.fieldnames), 4 if level == "t2r" else 3)
                    rows = list(reader)
                self.assertEqual(len(rows), 1)
                self.assertNotIn("group", rows[0])
                self.assertIn(f"{4 if level == 't2r' else 3} columns", stdout.getvalue())
                metric_values = [value for key, value in rows[0].items()
                                 if key not in extract.METADATA]
                self.assertEqual(metric_values, ["7", "5"] if level == "t2r" else ["3"])


if __name__ == "__main__":
    unittest.main()
