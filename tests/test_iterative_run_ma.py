"""Check cohort discovery, safe extractor calls and per-session CSV exports."""

import contextlib
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import nibabel as nib
import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "bin/3.2.1_DTIdata_extract/batch_DTIdata_extract.py"
spec = importlib.util.spec_from_file_location("iterative_run_ma", SCRIPT)
batch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(batch)
extract_spec = importlib.util.spec_from_file_location("batch_test_dti_extract", SCRIPT.with_name("DTIdata_extract.py"))
extract = importlib.util.module_from_spec(extract_spec)
extract_spec.loader.exec_module(extract)


class IterativeRunMATests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.root = self.folder / "proc folder"
        self.root.mkdir()
        self.translator = self.folder / "translations with spaces.txt"
        self.translator.write_text("22\tPTLp\n", encoding="utf-8")

    def session(self, subject="mouse A", session="PT 3", roi="scan_AnnoSplit_parental.nii.gz",
                metrics=("fa", "md", "rd", "ad"), extension=".nii.gz"):
        dwi = self.root / subject / session / "dwi"
        maps = dwi / "DSI_studio"
        maps.mkdir(parents=True)
        if roi is not None:
            (dwi / roi).touch()
        for metric in metrics:
            (maps / f"scan.fz.{metric}{extension}").touch()
        return dwi

    def run_mocked(self, roi_name="annosplit_parental", extra=(), side_effect=None):
        process_logs_before = {path: path.read_bytes() for path in self.root.rglob("*process.log")
                               if path.is_file()}
        logs_before = {path: path.read_bytes() for path in self.root.rglob("*analysis.log") if path.is_file()}
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            with patch.object(batch.subprocess, "run", return_value=SimpleNamespace(returncode=0),
                              side_effect=side_effect) as run:
                status = batch.main(["-i", str(self.root), "-r", roi_name, *map(str, extra)])
        logs_after = {path: path.read_bytes() for path in self.root.rglob("*analysis.log") if path.is_file()}
        self.assertEqual({path: path.read_bytes() for path in self.root.rglob("*process.log")
                          if path.is_file()}, process_logs_before)
        for path, contents in logs_after.items():
            previous = logs_before.get(path, b"")
            self.assertTrue(contents.startswith(previous))
            added = contents[len(previous):].decode("utf-8")
            if added:
                self.assertEqual(added.count("=== DTIdata_extract.py "), 1)
                self.assertNotIn("Image:", added)
                self.assertNotIn("Done:", added)
        self.assertTrue(all(line.startswith(("Summary:", "Failed sessions:", "Skipped sessions:", "  "))
                            for line in stdout.getvalue().splitlines()))
        return status, stdout.getvalue(), stderr.getvalue(), run

    def log_path(self, dwi):
        subject, session = dwi.parent.parent.name, dwi.parent.name
        name = f"{subject}_{session}_dwi_analysis.log" if subject.startswith("sub-") and session.startswith("ses-") else "analysis.log"
        return dwi / name

    def log_text(self, dwi=None):
        dwi = dwi or self.root / "mouse A" / "PT 3" / "dwi"
        return self.log_path(dwi).read_text(encoding="utf-8")

    def copied_driver_env(self):
        # Copied driver scripts must still be able to import the real common package.
        env = os.environ.copy()
        env["PYTHONPATH"] = str(SCRIPT.parents[1]) + os.pathsep + env.get("PYTHONPATH", "")
        return env

    def manifest(self, dwi):
        path = dwi / f".{dwi.parent.parent.name}_{dwi.parent.name}_dwi_aidamri_manifest.json"
        return json.loads(path.read_text())

    def write_export(self, command, rows):
        path = Path(batch.getOutfile(command[command.index("-r") + 1], command[command.index("-i") + 1],
                                     exclude_stroke="--exclude-stroke-mask" in command))
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["roi_id", "acronym", "median"])
            writer.writerows(rows)
        return SimpleNamespace(returncode=0)

    def read_table(self, name):
        with (self.root / name).open(encoding="utf-8", newline="") as stream:
            return list(csv.reader(stream))

    def test_stroke_mask_matches_case_insensitive_nifti_suffix_only_directly_in_dwi(self):
        dwi = self.session()
        (dwi / "nested").mkdir()
        for ignored in (dwi / "nested" / "nested_Stroke_mask.nii.gz",
                        dwi / "DSI_studio" / "nested_Stroke_mask.nii.gz",
                        dwi / "scan_Stroke_mask_anno.nii.gz", dwi / "scan_Stroke_mask.nii.txt"):
            ignored.touch()
        (dwi / "directory_Stroke_mask.nii.gz").mkdir()
        for name in ("Stroke_mask.nii", "scan_StRoKe_MaSk.NII.GZ", "scan_Stroke_mask.nii.gz"):
            mask = dwi / name
            mask.touch()
            try:
                with self.subTest(name=name):
                    status, _, stderr, run = self.run_mocked(extra=("--exclude-stroke-mask",))
                    self.assertEqual(status, 0, stderr)
                    self.assertEqual(stderr, "")
                    self.assertEqual(run.call_count, 4)
                    for call in run.call_args_list:
                        command = call.args[0]
                        self.assertEqual(command[command.index("--exclude-stroke-mask") + 1], str(mask))
                        self.assertTrue(Path(command[command.index("--exclude-stroke-mask") + 1]).is_absolute())
            finally:
                mask.unlink()

    def test_missing_stroke_mask_warns_once_and_forwards_no_mask_argument(self):
        dwi = self.session()
        (dwi / "DSI_studio" / "scan_Stroke_mask.nii.gz").touch()
        (dwi / "scan_Stroke_mask_anno.nii.gz").touch()
        status, stdout, stderr, run = self.run_mocked(extra=("--exclude-stroke-mask",))
        self.assertEqual(status, 0, stderr)
        self.assertEqual(run.call_count, 4)
        self.assertEqual(stderr.count("No stroke mask"), 1)
        self.assertIn("mouse A/PT 3/dwi", stderr)
        self.assertIn("Calculating medians without stroke exclusion", stderr)
        self.assertIn("4 succeeded, 0 skipped, 0 failed", stdout)
        for call in run.call_args_list:
            self.assertNotIn("--exclude-stroke-mask", call.args[0])

    def test_without_switch_does_not_search_or_forward_stroke_masks(self):
        dwi = self.session()
        (dwi / "a_Stroke_mask.nii").touch()
        (dwi / "b_Stroke_mask.nii.gz").touch()
        self.session("other subject", "other session")
        with patch.object(batch, "_find_rois", wraps=batch._find_rois) as find:
            status, _, stderr, run = self.run_mocked()
        self.assertEqual(status, 0, stderr)
        self.assertEqual(stderr, "")
        self.assertEqual(run.call_count, 8)
        self.assertTrue(all(call.args[1] == "annosplit_parental" for call in find.call_args_list))
        for call in run.call_args_list:
            self.assertNotIn("--exclude-stroke-mask", call.args[0])

    def test_multiple_stroke_masks_abort_entire_batch_before_outputs(self):
        first = self.session("sub-A", "ses-1", metrics=("fa",))
        (first / "scan_Stroke_mask.nii.gz").touch()
        ambiguous = self.session("sub-Z", "ses-2", metrics=("fa",))
        masks = [ambiguous / "a_Stroke_mask.nii", ambiguous / "b_STROKE_MASK.NII.GZ"]
        for mask in masks:
            mask.touch()
        output = self.folder / "new output folder" / "cohort"
        status, stdout, stderr, run = self.run_mocked(
            extra=("--exclude-stroke-mask", "-t", self.translator, "-o", output, "-m", "fa"),
        )
        self.assertEqual(status, 1)
        run.assert_not_called()
        self.assertIn("sub-Z/ses-2/dwi", stderr)
        self.assertIn("Multiple matching stroke masks", stderr)
        for mask in masks:
            self.assertIn(mask.name, stderr)
        self.assertIn("Aborting batch", stderr)
        self.assertIn("0 succeeded", stdout)
        self.assertIn("sub-Z / ses-2", stdout)
        self.assertFalse(output.parent.exists())
        self.assertEqual(list(self.root.rglob("*.csv")), [])
        self.assertEqual(list(self.root.rglob("*analysis.log")), [])

    def test_stroke_mask_search_error_aborts_before_any_extraction(self):
        self.session()
        with patch.object(batch, "_find_rois", side_effect=PermissionError("directory denied")):
            status, stdout, stderr, run = self.run_mocked(extra=("--exclude-stroke-mask",))
        self.assertEqual(status, 1)
        run.assert_not_called()
        self.assertIn("mouse A/PT 3/dwi", stderr)
        self.assertIn("Cannot search for stroke masks: directory denied", stderr)
        self.assertIn("0 succeeded", stdout)

    def test_stroke_mask_extraction_failure_does_not_retry_without_mask(self):
        dwi = self.session(metrics=("fa",))
        (dwi / "scan_Stroke_mask.nii.gz").touch()
        status, _, stderr, run = self.run_mocked(
            extra=("--exclude-stroke-mask", "-m", "fa"),
            side_effect=lambda *args, **kwargs: SimpleNamespace(returncode=1),
        )
        self.assertEqual(status, 1)
        self.assertEqual(run.call_count, 1)
        self.assertIn("--exclude-stroke-mask", run.call_args.args[0])
        self.assertIn("Extraction failed", stderr)
        self.assertNotIn("No stroke mask", stderr)

    def test_empty_medians_are_accepted_and_nonempty_values_must_be_finite_numbers(self):
        source = self.folder / "regions.csv"
        source.write_text("roi_id,acronym,median\n1,A,\n2,B,0.25\n", encoding="utf-8")
        self.assertEqual(batch._read_region_values(source), {"A": "", "B": "0.25"})
        for invalid in ("nan", "inf", "-inf", "invalid", " "):
            with self.subTest(invalid=invalid):
                source.write_text(f"roi_id,acronym,median\n1,A,{invalid}\n", encoding="utf-8")
                with self.assertRaises(ValueError):
                    batch._read_region_values(source)

    def test_end_to_end_stroke_mask_cohort_combines_masked_and_unmasked_sessions(self):
        sessions = [self.session(subject, "ses-1", metrics=("fa", "md"))
                    for subject in ("sub-A", "sub-B", "sub-C")]
        labels = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        image = np.array([np.nan, 0.1, 0.3, 0.7], dtype=np.float32).reshape(2, 2, 1)
        originals = {}
        stale_outputs = []
        for index, dwi in enumerate(sessions):
            roi = dwi / "scan_AnnoSplit_parental.nii.gz"
            nib.save(nib.Nifti1Image(labels, np.eye(4)), roi)
            originals[roi] = roi.read_bytes()
            if index < 2:
                mask = np.array([0, 0, 1, 1] if index == 0 else [1, 1, 1, 1], dtype=np.uint8)
                path = dwi / "scan_StRoKe_MaSk.nii.gz"
                nib.save(nib.Nifti1Image(mask.reshape(labels.shape), np.eye(4)), path)
                originals[path] = path.read_bytes()
            for metric in ("fa", "md"):
                path = dwi / "DSI_studio" / f"scan.fz.{metric}.nii.gz"
                nib.save(nib.Nifti1Image(image, np.eye(4)), path)
                originals[path] = path.read_bytes()
                stale = Path(batch.getOutfile(roi, path, exclude_stroke=index == 2))
                stale.write_bytes(b"Previous result from the other variant\n")
                stale_outputs.append(stale)
        result = subprocess.run(
            [sys.executable, "-B", str(SCRIPT), "-i", str(self.root), "-r", "annosplit_parental",
             "-t", str(self.translator), "-o", str(self.root / "cohort.csv"), "--exclude-stroke-mask"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stderr.count("No stroke mask"), 1)
        self.assertIn("sub-C/ses-1/dwi", result.stderr)
        for metric in ("fa", "md"):
            for subject, roi_ids in (("sub-A", (2022,)), ("sub-B", (22, 2022))):
                for roi_id in roi_ids:
                    warning = (f"Warning: [subject={subject}, session=ses-1, metric={metric.upper()}] "
                               f"All voxels in ROI ID {roi_id} excluded by stroke mask. Leaving median empty.")
                    self.assertIn(warning, result.stderr)
                    self.assertIn(warning, self.log_text(self.root / subject / "ses-1" / "dwi"))
            table = self.read_table(f"cohort_{metric}.csv")
            self.assertEqual(table[0], ["subject", "session", "L_PTLp", "R_PTLp"])
            self.assertEqual(len(table), 4)
            for index, (row, dwi) in enumerate(zip(table[1:], sessions)):
                self.assertEqual(row[:2], [dwi.parent.parent.name, "ses-1"])
                if index == 0:
                    self.assertAlmostEqual(float(row[2]), 0.1)
                    self.assertEqual(row[3], "")
                elif index == 1:
                    self.assertEqual(row[2:], ["", ""])
                else:
                    self.assertAlmostEqual(float(row[2]), 0.2)
                    self.assertAlmostEqual(float(row[3]), 0.7)
                source = Path(batch.getOutfile(dwi / "scan_AnnoSplit_parental.nii.gz",
                                              dwi / "DSI_studio" / f"scan.fz.{metric}.nii.gz",
                                              exclude_stroke=index < 2))
                with source.open(newline="") as stream:
                    values = {entry["acronym"]: entry["median"] for entry in csv.DictReader(stream)}
                self.assertEqual(row[2:], [values[acronym] for acronym in table[0][2:]])
                self.assertIn(source.relative_to(dwi).as_posix(),
                              self.manifest(dwi)["stages"]["extraction"]["created_files"])
            self.assertFalse((self.root / f"cohort_{metric}_exclude-stroke.csv").exists())
        for path, contents in originals.items():
            self.assertEqual(path.read_bytes(), contents)
        for stale in stale_outputs:
            self.assertEqual(stale.read_bytes(), b"Previous result from the other variant\n")

    def test_groups_are_added_after_session_for_all_metrics_and_sessions(self):
        for subject, session in (("sub-A", "ses-2"), ("sub-A", "ses-1"),
                                 ("sub-B", "ses-1"), ("sub-unknown", "ses-1")):
            self.session(subject, session, metrics=("fa", "md"))
        group_file = self.folder / "groups with spaces.csv"
        group_file.write_text(' tDCS ,"Sham, control"\n sub-A , sub-B \n,sub-outside\n', encoding="utf-8")

        status, _, stderr, run = self.run_mocked(
            extra=("-t", self.translator, "-o", str(self.root / "cohort"), "-g", group_file),
            side_effect=lambda command, **kwargs: self.write_export(
                command, [(1, "Z", "0.1234567890123456789"), (2, "A", "0.5")],
            ),
        )
        self.assertEqual(status, 0, stderr)
        self.assertEqual(run.call_count, 8)
        expected = [
            ["subject", "session", "group", "Z", "A"],
            ["sub-A", "ses-1", "tDCS", "0.1234567890123456789", "0.5"],
            ["sub-A", "ses-2", "tDCS", "0.1234567890123456789", "0.5"],
            ["sub-B", "ses-1", "Sham, control", "0.1234567890123456789", "0.5"],
            ["sub-unknown", "ses-1", "", "0.1234567890123456789", "0.5"],
        ]
        for metric in ("fa", "md"):
            self.assertEqual(self.read_table(f"cohort_{metric}.csv"), expected)
        for call in run.call_args_list:
            self.assertNotIn("-g", call.args[0])
            self.assertNotIn("--group", call.args[0])

    def test_group_csv_accepts_bom_and_common_delimiters_with_unequal_group_sizes(self):
        self.session("sub-A", "ses-1", metrics=("fa",))
        self.session("sub-B", "ses-1", metrics=("fa",))
        group_file = self.folder / "groups.csv"
        for delimiter in (",", ";", "\t"):
            for bom in (False, True):
                with self.subTest(delimiter=delimiter, bom=bom):
                    with group_file.open("w", encoding="utf-8-sig" if bom else "utf-8", newline="") as stream:
                        writer = csv.writer(stream, delimiter=delimiter)
                        writer.writerows([
                            ["tDCS", "Sham"], ["sub-A", "sub-B"], ["sub-outside"], [], ["", "sub-other"],
                        ])
                    status, _, stderr, _ = self.run_mocked(
                        extra=("-t", self.translator, "-o", str(self.root / "cohort"), "--group", group_file, "-m", "fa"),
                        side_effect=lambda command, **kwargs: self.write_export(command, [(1, "A", "0.5")]),
                    )
                    self.assertEqual(status, 0, stderr)
                    self.assertEqual(self.read_table("cohort_fa.csv"), [
                        ["subject", "session", "group", "A"],
                        ["sub-A", "ses-1", "tDCS", "0.5"], ["sub-B", "ses-1", "Sham", "0.5"],
                    ])

    def test_invalid_group_csv_fails_before_extraction_and_preserves_previous_output(self):
        self.session("sub-A", "ses-1", metrics=("fa",))
        group_file = self.folder / "groups.csv"
        output = self.root / "cohort_fa.csv"
        output.write_bytes(b"Previous cohort\n")
        invalid = (
            b"", b"\xef\xbb\xbf", b",Sham\nsub-A,sub-B\n", b"tDCS, \nsub-A,sub-B\n",
            b"tDCS,tDCS\nsub-A,sub-B\n", b"tDCS, tDCS \nsub-A,sub-B\n",
            b"tDCS,Sham\nsub-A,sub-A\n", b"tDCS,Sham\nsub-A,sub-B\nsub-A,\n",
            b"tDCS,Sham\nsub-A,sub-B,sub-C\n", b'tDCS,Sham\n"unterminated,sub-B\n',
            b"tDCS,Sham\n\xff,sub-B\n",
        )
        for contents in invalid:
            with self.subTest(contents=contents):
                group_file.write_bytes(contents)
                status, stdout, _, run = self.run_mocked(
                    extra=("-t", self.translator, "-o", str(self.root / "cohort"), "-g", group_file),
                )
                self.assertEqual(status, 1)
                self.assertIn("Cannot read group CSV", stdout)
                run.assert_not_called()
                self.assertEqual(output.read_bytes(), b"Previous cohort\n")
        self.assertEqual(list(self.root.rglob("*analysis.log")), [])

    def test_group_option_requires_output_and_an_existing_csv(self):
        self.session(metrics=("fa",))
        group_file = self.folder / "groups.csv"
        group_file.write_text("tDCS,Sham\nsub-A,sub-B\n")
        status, stdout, _, run = self.run_mocked(extra=("-g", group_file))
        self.assertEqual(status, 1)
        self.assertIn("--group requires -o/--output", stdout)
        run.assert_not_called()
        for path in (self.folder / "absent.csv", self.folder):
            with self.subTest(path=path):
                status, stdout, _, run = self.run_mocked(
                    extra=("-t", self.translator, "-o", str(self.root / "cohort"), "-g", path),
                )
                self.assertEqual(status, 1)
                self.assertIn("Required file does not exist", stdout)
                run.assert_not_called()
        self.assertEqual(list(self.root.rglob("*analysis.log")), [])

    def test_group_column_remains_present_when_no_subject_has_a_matching_group(self):
        self.session("sub-NR1037", "ses-1", metrics=("fa",))
        group_file = self.folder / "groups.csv"
        for contents in ("tDCS,Sham\n", "tDCS,Sham\n,\n", "tDCS,Sham\nNR1037,SUB-NR1037\n"):
            with self.subTest(contents=contents):
                group_file.write_text(contents)
                status, _, stderr, _ = self.run_mocked(
                    extra=("-t", self.translator, "-o", str(self.root / "cohort"), "-g", group_file, "-m", "fa"),
                    side_effect=lambda command, **kwargs: self.write_export(command, [(1, "A", "0.5")]),
                )
                self.assertEqual(status, 0, stderr)
                self.assertEqual(stderr, "")
                self.assertEqual(self.read_table("cohort_fa.csv"), [
                    ["subject", "session", "group", "A"], ["sub-NR1037", "ses-1", "", "0.5"],
                ])

    def test_group_acronym_conflicts_only_when_group_metadata_is_enabled(self):
        self.session("sub-A", "ses-1", metrics=("fa",))
        group_file = self.folder / "groups.csv"
        group_file.write_text("tDCS\nsub-A\n")
        export = lambda command, **kwargs: self.write_export(command, [(1, "group", "0.5")])
        arguments = ("-t", self.translator, "-o", str(self.root / "cohort"), "-m", "fa")
        status, _, stderr, _ = self.run_mocked(extra=(*arguments, "-g", group_file), side_effect=export)
        self.assertEqual(status, 1)
        self.assertIn("reserved acronym 'group'", stderr)
        self.assertFalse((self.root / "cohort_fa.csv").exists())
        status, _, stderr, _ = self.run_mocked(extra=arguments, side_effect=export)
        self.assertEqual(status, 0, stderr)
        self.assertEqual(self.read_table("cohort_fa.csv"), [
            ["subject", "session", "group"], ["sub-A", "ses-1", "0.5"],
        ])

    def test_end_to_end_group_column_annotates_real_exports_and_leaves_unknown_subject_blank(self):
        sessions = [self.session("sub-A", "ses-1", metrics=("fa",)),
                    self.session("sub-A", "ses-2", metrics=("fa",)),
                    self.session("sub-unknown", "ses-1", metrics=("fa",))]
        group_file = self.folder / "groups.csv"
        group_file.write_text("tDCS,Sham\nsub-A,sub-outside\n")
        labels = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        values = np.array([np.nan, 0.1, 0.3, 0.7], dtype=np.float32).reshape(2, 2, 1)
        for dwi in sessions:
            nib.save(nib.Nifti1Image(labels, np.eye(4)), dwi / "scan_AnnoSplit_parental.nii.gz")
            nib.save(nib.Nifti1Image(values, np.eye(4)), dwi / "DSI_studio" / "scan.fz.fa.nii.gz")
        result = subprocess.run(
            [sys.executable, "-B", str(SCRIPT), "-i", str(self.root), "-r", "annosplit_parental",
             "-t", str(self.translator), "-o", str(self.root / "cohort.csv"), "--group", str(group_file), "-m", "fa"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stderr, "")
        table = self.read_table("cohort_fa.csv")
        self.assertEqual(table[0], ["subject", "session", "group", "L_PTLp", "R_PTLp"])
        for row, dwi in zip(table[1:], sessions):
            self.assertEqual(row[:3], [dwi.parent.parent.name, dwi.parent.name,
                                      "tDCS" if dwi.parent.parent.name == "sub-A" else ""])
            with (dwi / "DSI_studio" / "scan.fz.fa_AnnoSplit_parental.csv").open(newline="") as stream:
                source = {entry["acronym"]: entry["median"] for entry in csv.DictReader(stream)}
            self.assertEqual(row[3:], [source[acronym] for acronym in table[0][3:]])

    def test_cohort_output_joins_acronyms_and_preserves_median_text(self):
        self.session("sub-B", "ses-1", metrics=("fa", "md"))
        self.session("sub-A", "ses-2", metrics=("fa", "md"))
        self.session("sub-A", "ses-1", metrics=("fa", "md"))
        special = 'Region, "A"'

        def export(command, **kwargs):
            image = Path(command[command.index("-i") + 1])
            if ".md." in image.name:
                rows = [(9, "MD region", "-0.000000123456789012345")]
            elif kwargs["cwd"].parent.parent.name == "sub-B":
                rows = [(3, "B", "1.23000000000000001e-8"), (4, special, "0.000"), (1, "Z", "0.25")]
            else:
                rows = [(1, "Z", "0.1234567890123456789"), (2, "A", "0.5")]
            return self.write_export(command, rows)

        status, stdout, stderr, run = self.run_mocked(
            extra=("-t", self.translator, "-o", str(self.root / "cohort.csv")), side_effect=export,
        )
        self.assertEqual(status, 0, stderr)
        self.assertEqual(run.call_count, 6)
        self.assertIn(f"Cohort CSVs: {self.root / 'cohort_fa.csv'}, {self.root / 'cohort_md.csv'}.", stdout)
        self.assertEqual(self.read_table("cohort_fa.csv"), [
            ["subject", "session", "Z", "A", "B", special],
            ["sub-A", "ses-1", "0.1234567890123456789", "0.5", "", ""],
            ["sub-A", "ses-2", "0.1234567890123456789", "0.5", "", ""],
            ["sub-B", "ses-1", "0.25", "", "1.23000000000000001e-8", "0.000"],
        ])
        self.assertEqual(self.read_table("cohort_md.csv"), [
            ["subject", "session", "MD region"],
            ["sub-A", "ses-1", "-0.000000123456789012345"],
            ["sub-A", "ses-2", "-0.000000123456789012345"],
            ["sub-B", "ses-1", "-0.000000123456789012345"],
        ])
        self.assertIn("No valid RD rows", stderr)
        self.assertIn("No valid AD rows", stderr)
        self.assertFalse((self.root / "cohort.csv").exists())
        self.assertEqual(len(list(self.root.glob("*.csv"))), 2)

    def test_output_requires_translation_and_a_base_filename_before_any_extraction(self):
        self.session(metrics=("fa",))
        status, stdout, _, run = self.run_mocked(extra=("-o", str(self.root / "cohort.csv")))
        self.assertEqual(status, 1)
        self.assertIn("requires -t/--translatorTXT", stdout)
        run.assert_not_called()
        for name in ("", " ", ".", "..", ".csv", "nested/.csv", "nested/.", "nested/..",
                     "nested/", "bad\0name", str(self.root)):
            with self.subTest(name=name):
                status, stdout, _, run = self.run_mocked(extra=("-t", self.translator, "-o", name))
                self.assertEqual(status, 1)
                self.assertIn("Error:", stdout)
                run.assert_not_called()
        self.assertEqual(list(self.root.rglob("*analysis.log")), [])

    def test_output_paths_create_directories_and_append_each_metric(self):
        self.session()
        for name in ("cohort", "cohort.csv"):
            output_base = self.folder / "new output" / name / name
            with self.subTest(name=name):
                status, stdout, stderr, run = self.run_mocked(
                    extra=("-t", self.translator, "-o", output_base),
                    side_effect=lambda command, **kwargs: self.write_export(command, [(1, "A", "0.25")]),
                )
                self.assertEqual(status, 0, stderr)
                self.assertEqual(run.call_count, 4)
                for metric in batch.METRICS:
                    output = output_base.parent / f"cohort_{metric}.csv"
                    self.assertEqual(self.read_table(output), [
                        ["subject", "session", "A"], ["mouse A", "PT 3", "0.25"],
                    ])
                    self.assertIn(str(output), stdout)
                self.assertEqual(len(list(output_base.parent.glob("*.csv"))), 4)
        self.assertEqual(list(self.root.glob("*.csv")), [])

    def test_relative_and_home_output_paths_use_invocation_directory(self):
        self.session(metrics=("fa",))
        working_directory = self.folder / "invocation directory"
        working_directory.mkdir()
        previous_directory = Path.cwd()
        try:
            os.chdir(working_directory)
            with patch.dict(os.environ, {"HOME": str(self.folder)}):
                for name, expected in (
                    ("cohort", working_directory / "cohort_fa.csv"),
                    ("nested/cohort.csv", working_directory / "nested/cohort_fa.csv"),
                    ("../outside/cohort", self.folder / "outside/cohort_fa.csv"),
                    ("~/home output/cohort", self.folder / "home output/cohort_fa.csv"),
                ):
                    with self.subTest(name=name):
                        status, stdout, stderr, _ = self.run_mocked(
                            extra=("-t", self.translator, "-o", name, "-m", "fa"),
                            side_effect=lambda command, **kwargs: self.write_export(command, [(1, "A", "0.25")]),
                        )
                        self.assertEqual(status, 0, stderr)
                        self.assertTrue(expected.is_file())
                        self.assertIn(str(expected), stdout)
        finally:
            os.chdir(previous_directory)
        self.assertEqual(list(self.root.glob("*.csv")), [])

    def test_invalid_output_destinations_fail_before_extraction(self):
        self.session(metrics=("fa",))
        blocked_parent = self.folder / "blocked"
        blocked_parent.write_bytes(b"Existing file\n")
        directory_target = self.folder / "cohort_fa.csv"
        directory_target.mkdir()
        for name in (blocked_parent / "cohort", self.folder / "cohort"):
            with self.subTest(name=name):
                status, stdout, _, run = self.run_mocked(
                    extra=("-t", self.translator, "-o", name, "-m", "fa"),
                )
                self.assertEqual(status, 1)
                self.assertIn("Cannot prepare cohort output", stdout)
                run.assert_not_called()
        self.assertEqual(blocked_parent.read_bytes(), b"Existing file\n")
        self.assertTrue(directory_target.is_dir())

    def test_output_directory_creation_failure_aborts_before_extraction(self):
        self.session(metrics=("fa",))
        with patch.object(batch.Path, "mkdir", side_effect=PermissionError("mkdir denied")):
            status, stdout, _, run = self.run_mocked(
                extra=("-t", self.translator, "-o", self.folder / "new/cohort", "-m", "fa"),
            )
        self.assertEqual(status, 1)
        self.assertIn("mkdir denied", stdout)
        run.assert_not_called()
        self.assertFalse((self.folder / "new").exists())

    def test_output_accepts_long_option_spaces_and_missing_or_uppercase_extension(self):
        self.session(metrics=("fa", "md"))
        for name in ("cohort result", "cohort result.CSV"):
            with self.subTest(name=name):
                status, _, stderr, run = self.run_mocked(
                    extra=("-t", self.translator, "--output", self.root / name, "-m", "md"),
                    side_effect=lambda command, **kwargs: self.write_export(command, [(1, "A", "0.25")]),
                )
                self.assertEqual(status, 0, stderr)
                self.assertEqual(run.call_count, 1)
                self.assertEqual(self.read_table("cohort result_md.csv"), [
                    ["subject", "session", "A"], ["mouse A", "PT 3", "0.25"],
                ])
        self.assertEqual([p.name for p in self.root.glob("*.csv")], ["cohort result_md.csv"])

    def test_cohort_excludes_stale_failed_skipped_and_unselected_atlas_results(self):
        good = self.session("sub-good", "ses-1", metrics=("fa", "md"))
        bad = self.session("sub-bad", "ses-1", metrics=("fa", "md"))
        skipped = self.session("sub-skipped", "ses-1", roi=None, metrics=("fa",))
        stale = bad / "DSI_studio" / "scan.fz.fa_AnnoSplit_parental.csv"
        stale.write_text("roi_id,acronym,median\n1,STALE,99\n")
        (skipped / "DSI_studio" / stale.name).write_text(stale.read_text())
        (good / "DSI_studio" / "scan.fz.fa_Anno.csv").write_text(stale.read_text())

        def export(command, **kwargs):
            image = Path(command[command.index("-i") + 1])
            if kwargs["cwd"] == bad and ".fa." in image.name:
                return SimpleNamespace(returncode=1)
            return self.write_export(command, [(1, "CURRENT", "0.5")])

        status, stdout, _, _ = self.run_mocked(
            extra=("-t", self.translator, "-o", str(self.root / "cohort")), side_effect=export,
        )
        self.assertEqual(status, 1)
        self.assertIn("sub-bad / ses-1", stdout)
        self.assertEqual(self.read_table("cohort_fa.csv"), [
            ["subject", "session", "CURRENT"], ["sub-good", "ses-1", "0.5"],
        ])
        self.assertEqual(self.read_table("cohort_md.csv"), [
            ["subject", "session", "CURRENT"],
            ["sub-bad", "ses-1", "0.5"], ["sub-good", "ses-1", "0.5"],
        ])
        self.assertIn("STALE", stale.read_text())

    def test_ambiguous_fa_and_fa0_are_extracted_but_omitted_from_cohort(self):
        ambiguous = self.session("sub-A", "ses-1", metrics=("fa", "fa0", "md"))
        self.session("sub-B", "ses-1", metrics=("fa", "md"))
        export = lambda command, **kwargs: self.write_export(command, [(1, "A", "0.1")])
        status, stdout, stderr, run = self.run_mocked(
            extra=("-t", self.translator, "-o", str(self.root / "cohort")), side_effect=export,
        )
        self.assertEqual(status, 1)
        self.assertEqual(run.call_count, 5)
        self.assertIn("Multiple FA images", stderr)
        self.assertIn("sub-A / ses-1", stdout)
        self.assertEqual(self.read_table("cohort_fa.csv"), [
            ["subject", "session", "A"], ["sub-B", "ses-1", "0.1"],
        ])
        self.assertEqual(len(self.read_table("cohort_md.csv")), 3)
        self.assertEqual(len(list((ambiguous / "DSI_studio").glob("*.csv"))), 3)
        # An unselected FA0 must not make an explicit -m fa ambiguous.
        status, _, stderr, _ = self.run_mocked(
            extra=("-t", self.translator, "-o", str(self.root / "cohort"), "-m", "fa"), side_effect=export,
        )
        self.assertEqual(status, 0, stderr)
        self.assertEqual(len(self.read_table("cohort_fa.csv")), 3)

    def test_bad_export_csvs_fail_collection_without_replacing_previous_cohort(self):
        dwi = self.session(metrics=("fa",))
        source = dwi / "DSI_studio" / "scan.fz.fa_AnnoSplit_parental.csv"
        previous = self.root / "cohort_fa.csv"
        previous.write_bytes(b"Previous cohort\n")
        invalid = (
            b"wrong,header\n",
            b"roi_id,acronym,median\n1,A\n",
            b"roi_id,acronym,median\ninvalid,A,1\n",
            b"roi_id,acronym,median\n1,A,1\n2,A,2\n",
            b"roi_id,acronym,median\n1,A,1\n1,B,2\n",
            b"roi_id,acronym,median\n1,,1\n",
            b"roi_id,acronym,median\n1,subject,1\n",
            b"roi_id,acronym,median\n1,session,1\n",
            b"roi_id,acronym,median\n1,A,nan\n",
            b"roi_id,acronym,median\n1,A,inf\n",
            b"roi_id,acronym,median\n1,A,invalid\n",
            b'roi_id,acronym,median\n1,"unterminated\n',
            b"roi_id,acronym,median\n1,\xff,1\n",
            None,
        )
        for contents in invalid:
            def export(command, **kwargs):
                source.unlink(missing_ok=True)
                if contents is not None:
                    source.write_bytes(contents)
                return SimpleNamespace(returncode=0)

            with self.subTest(contents=contents):
                status, stdout, stderr, _ = self.run_mocked(
                    extra=("-t", self.translator, "-o", str(self.root / "cohort"), "-m", "fa"), side_effect=export,
                )
                self.assertEqual(status, 1)
                self.assertIn("Cannot collect", stderr)
                self.assertIn("mouse A / PT 3", stdout)
                self.assertEqual(previous.read_bytes(), b"Previous cohort\n")

    def test_empty_export_warns_and_leaves_existing_cohort_unchanged(self):
        self.session(metrics=("fa",))
        previous = self.root / "cohort_fa.csv"
        previous.write_bytes(b"Previous cohort\n")
        status, _, stderr, _ = self.run_mocked(
            extra=("-t", self.translator, "-o", str(self.root / "cohort"), "-m", "fa"),
            side_effect=lambda command, **kwargs: self.write_export(command, []),
        )
        self.assertEqual(status, 0, stderr)
        self.assertIn("No region values", stderr)
        self.assertIn("No valid FA rows", stderr)
        self.assertEqual(previous.read_bytes(), b"Previous cohort\n")

    def test_failed_atomic_replace_preserves_existing_cohort_and_removes_temporary_file(self):
        self.session(metrics=("fa",))
        previous = self.root / "cohort_fa.csv"
        previous.write_bytes(b"Previous cohort\n")
        with patch.object(batch.os, "replace", side_effect=PermissionError("replace denied")):
            status, stdout, stderr, _ = self.run_mocked(
                extra=("-t", self.translator, "-o", str(self.root / "cohort"), "-m", "fa"),
                side_effect=lambda command, **kwargs: self.write_export(command, [(1, "A", "0.5")]),
            )
        self.assertEqual(status, 1)
        self.assertIn("1 succeeded, 0 skipped, 1 failed", stdout)
        self.assertIn("Cannot write cohort CSV", stderr)
        self.assertIn("replace denied", stderr)
        self.assertEqual(previous.read_bytes(), b"Previous cohort\n")
        self.assertEqual(list(self.root.glob(".*.tmp")), [])

    def test_cohort_permissions_match_normal_file_creation_under_umask(self):
        output = self.root / "cohort_fa.csv"
        rows = {("sub-A", "ses-1"): {"A": "0.5"}}
        for mask in (0o022, 0o027, 0o077):
            with self.subTest(umask=oct(mask)):
                reference = self.root / f"reference_{mask}.csv"
                original = os.umask(mask)
                try:
                    with reference.open("x", encoding="utf-8") as stream:
                        stream.write("reference\n")
                    batch._write_cohort_csv(output, rows)
                finally:
                    os.umask(original)
                self.assertEqual(stat.S_IMODE(output.stat().st_mode), stat.S_IMODE(reference.stat().st_mode))
                self.assertEqual(self.read_table(output.name), [
                    ["subject", "session", "A"], ["sub-A", "ses-1", "0.5"],
                ])
        self.assertEqual(list(self.root.glob(".*.tmp")), [])

    def test_cohort_temporary_name_collision_preserves_existing_file(self):
        output = self.root / "cohort_fa.csv"
        collision = self.root / f".{output.name}.collision.tmp"
        collision.write_bytes(b"Existing unrelated file\n")
        identifiers = [SimpleNamespace(hex="collision"), SimpleNamespace(hex="new")]
        with patch.object(batch.uuid, "uuid4", side_effect=identifiers):
            batch._write_cohort_csv(output, {("sub-A", "ses-1"): {"A": "0.5"}})
        self.assertEqual(collision.read_bytes(), b"Existing unrelated file\n")
        self.assertEqual(self.read_table(output.name), [
            ["subject", "session", "A"], ["sub-A", "ses-1", "0.5"],
        ])
        self.assertFalse((self.root / f".{output.name}.new.tmp").exists())

    def test_cohort_is_rebuilt_from_current_rows_on_repeated_runs(self):
        self.session("sub-A", "ses-1", metrics=("fa",))
        later = self.session("sub-B", "ses-1", metrics=("fa",))
        arguments = ("-t", self.translator, "-o", str(self.root / "cohort"), "-m", "fa")
        status, _, stderr, _ = self.run_mocked(
            extra=arguments,
            side_effect=lambda command, **kwargs: self.write_export(command, [(1, "OLD", "0.5")]),
        )
        self.assertEqual(status, 0, stderr)
        self.assertEqual(len(self.read_table("cohort_fa.csv")), 3)

        def updated(command, **kwargs):
            if kwargs["cwd"] == later:
                return SimpleNamespace(returncode=2)
            return self.write_export(command, [(2, "NEW", "0.75")])

        status, _, _, _ = self.run_mocked(extra=arguments, side_effect=updated)
        self.assertEqual(status, 1)
        self.assertEqual(self.read_table("cohort_fa.csv"), [
            ["subject", "session", "NEW"], ["sub-A", "ses-1", "0.75"],
        ])

    def test_end_to_end_cohort_values_equal_the_individual_csv_values(self):
        sessions = [self.session("sub-B", "ses-1", metrics=("fa", "md")),
                    self.session("sub-A", "ses-2", metrics=("fa", "md")),
                    self.session("sub-C", "ses-1", metrics=("fa", "md"))]
        self.translator.write_text("31\tACA\n22\tPTLp\n", encoding="utf-8")
        process_logs = []
        (self.root / "AIDAmri_git_information.txt").write_text("Cohort provenance\n")
        for index, dwi in enumerate(sessions):
            process_log = dwi / f"{dwi.parent.parent.name}_{dwi.parent.name}_dwi_process.log"
            process_log.write_bytes(b"Existing pipeline log\n")
            process_logs.append(process_log)
            labels = np.array([0, 22, 22, 2022 if index == 0 else 0], dtype=np.int16).reshape(2, 2, 1)
            if index == 2:
                labels.fill(0)
            nib.save(nib.Nifti1Image(labels, np.eye(4)), dwi / "scan_AnnoSplit_parental.nii.gz")
            for metric in ("fa", "md"):
                values = np.array([np.nan, 0.1 + index, 0.3 + index, 0.7], dtype=np.float32).reshape(2, 2, 1)
                nib.save(nib.Nifti1Image(values, np.eye(4)), dwi / "DSI_studio" / f"scan.fz.{metric}.nii.gz")
        result = subprocess.run(
            [sys.executable, "-B", str(SCRIPT), "-i", str(self.root), "-r", "annosplit_parental",
             "-t", str(self.translator), "-o", str(self.root / "cohort.csv")], capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for metric in ("fa", "md"):
            table = self.read_table(f"cohort_{metric}.csv")
            self.assertEqual(table[0], ["subject", "session", "L_ACA", "L_PTLp", "R_ACA", "R_PTLp"])
            self.assertEqual(len(table), 4)
            self.assertEqual(table[-1], ["sub-C", "ses-1", "", "", "", ""])
            for row, dwi in zip(table[1:], sorted(sessions)):
                self.assertEqual(row[2], "")
                self.assertEqual(row[4], "")
                self.assertEqual(row[:2], [dwi.parent.parent.name, dwi.parent.name])
                with (dwi / "DSI_studio" / f"scan.fz.{metric}_AnnoSplit_parental.csv").open(newline="") as stream:
                    source = {entry["acronym"]: entry["median"] for entry in csv.DictReader(stream)}
                self.assertEqual(list(source), table[0][2:])
                self.assertEqual(row[2:], [source[acronym] for acronym in table[0][2:]])
                stage = self.manifest(dwi)["stages"]["extraction"]
                self.assertIn(f"DSI_studio/scan.fz.{metric}_AnnoSplit_parental.csv", stage["created_files"])
        for dwi in sessions:
            self.assertEqual(self.log_text(dwi).count("=== DTIdata_extract.py "), 1)
            self.assertEqual(self.log_text(dwi).count("Cohort provenance"), 1)
            self.assertNotIn("Existing pipeline log", self.log_text(dwi))
        for process_log in process_logs:
            self.assertEqual(process_log.read_bytes(), b"Existing pipeline log\n")

    def test_discovers_sorted_subject_sessions_without_prefix_restrictions(self):
        last = self.session("subject-z", "session-z")
        first = self.session("animal-a", "time-a")
        middle = self.session("animal-a", "time-b")
        (self.root / "animal-a" / "time-a" / "anat").mkdir()
        (self.root / "shallow" / "dwi").mkdir(parents=True)
        (self.root / "too" / "deep" / "session" / "dwi").mkdir(parents=True)
        (self.root / "not" / "a-dir").mkdir(parents=True)
        (self.root / "not" / "a-dir" / "dwi").touch()
        self.assertEqual(batch.findData(self.root), [first, middle, last])

    def test_exact_case_insensitive_roi_suffix_and_both_nifti_extensions(self):
        dwi = self.session(roi="scan_aNnOsPlIt.NII.GZ", metrics=("fa",))
        parental = dwi / "scan_AnnoSplit_parental.nii"
        parental.touch()
        (dwi / "scan_AnnoSplit_extra.nii.gz").touch()
        (dwi / "scan_AnnoSplit.nii.txt").touch()
        (dwi / "nested_AnnoSplit.nii").mkdir()
        for name, expected in (("ANNOSPLIT", "scan_aNnOsPlIt.NII.GZ"),
                               ("annosplit_parental", parental.name)):
            with self.subTest(name=name):
                status, _, _, run = self.run_mocked(roi_name=name)
                self.assertEqual(status, 0)
                self.assertEqual(run.call_count, 1)
                command = run.call_args.args[0]
                self.assertEqual(command[command.index("-r") + 1], str(dwi / expected))

    def test_invokes_all_metrics_with_absolute_paths_and_current_interpreter(self):
        dwi = self.session(metrics=("fa", "fa0", "md", "rd", "ad"))
        (dwi / "DSI_studio" / "other.ad.nii").touch()
        (dwi / "DSI_studio" / "scan.qa.nii.gz").touch()
        (dwi / "DSI_studio" / "scan.fa.nii.gz.txt").touch()
        (dwi / "DSI_studio" / "directory.fa.nii.gz").mkdir()
        status, stdout, stderr, run = self.run_mocked(extra=("-t", self.translator))
        self.assertEqual(status, 0, stderr)
        self.assertEqual(stderr, "")
        self.assertEqual(run.call_count, 6)
        self.assertIn("Summary: 6 succeeded, 0 skipped, 0 failed.", stdout)
        self.assertNotIn("mouse A/PT 3/dwi", stdout)
        for call in run.call_args_list:
            command = call.args[0]
            self.assertEqual(command[0], sys.executable)
            self.assertEqual(Path(command[1]), SCRIPT.with_name("DTIdata_extract.py"))
            for option in ("-i", "-r", "-t"):
                self.assertTrue(Path(command[command.index(option) + 1]).is_absolute())
            self.assertEqual(command[-2:], ["-t", str(self.translator)])
            self.assertEqual(call.kwargs["cwd"], dwi)
            self.assertFalse(call.kwargs["check"])
            self.assertEqual(call.kwargs["stdout"], subprocess.DEVNULL)
            self.assertIsNone(call.kwargs["stderr"])

    def test_no_translation_option_is_forwarded_when_unspecified(self):
        self.session(metrics=("fa0",), extension=".nii")
        status, _, stderr, run = self.run_mocked()
        self.assertEqual(status, 0)
        self.assertNotIn("-t", run.call_args.args[0])
        self.assertIn("No MD map found", stderr)
        self.assertIn("No RD map found", stderr)
        self.assertIn("No AD map found", stderr)

    def test_single_metric_selection_matches_exact_suffix_and_both_extensions(self):
        dwi = self.session(metrics=(*batch.METRICS, "fa0"))
        for metric in batch.METRICS:
            (dwi / "DSI_studio" / f"other.{metric}.nii").touch()
        (dwi / "DSI_studio" / "scan.gfa.nii.gz").touch()
        (dwi / "DSI_studio" / "scan.fa.nii.gz.txt").touch()
        for metric in batch.METRICS:
            for option in ("-m", "--metric"):
                with self.subTest(metric=metric, option=option):
                    status, stdout, stderr, run = self.run_mocked(extra=(option, metric.upper()))
                    self.assertEqual(status, 0, stderr)
                    self.assertEqual(stderr, "")
                    self.assertEqual(run.call_count, 2)
                    names = {Path(call.args[0][3]).name for call in run.call_args_list}
                    self.assertEqual(names, {f"other.{metric}.nii", f"scan.fz.{metric}.nii.gz"})
                    self.assertIn("2 succeeded, 0 skipped, 0 failed", stdout)

    def test_explicit_all_matches_default_including_legacy_fa0(self):
        self.session(metrics=(*batch.METRICS, "fa0"))
        default_status, _, _, default_run = self.run_mocked()
        for option, value in (("-m", "all"), ("--metric", "ALL")):
            with self.subTest(option=option, value=value):
                status, stdout, stderr, run = self.run_mocked(extra=(option, value))
                self.assertEqual(status, default_status)
                self.assertEqual(stderr, "")
                self.assertEqual([call.args for call in run.call_args_list],
                                 [call.args for call in default_run.call_args_list])
                self.assertEqual(run.call_count, 5)
                self.assertIn("5 succeeded, 0 skipped, 0 failed", stdout)

    def test_missing_selected_metric_only_warns_for_that_metric(self):
        self.session(metrics=("fa",))
        status, stdout, stderr, run = self.run_mocked(extra=("-m", "md"))
        self.assertEqual(status, 1)
        run.assert_not_called()
        self.assertEqual(stderr.count("Warning:"), 1)
        self.assertIn("No MD map found", stderr)
        self.assertIn("0 succeeded, 1 skipped, 0 failed", stdout)
        self.assertIn("mouse A / PT 3 (fully skipped: MD missing)", stdout)

    def test_unselected_metric_collisions_do_not_affect_selected_metric(self):
        dwi = self.session(metrics=("fa", "ad"))
        (dwi / "DSI_studio" / "scan.fz.fa.nii").touch()
        status, stdout, stderr, run = self.run_mocked(extra=("-m", "ad"))
        self.assertEqual(status, 0)
        self.assertEqual(stderr, "")
        self.assertEqual(run.call_count, 1)
        self.assertTrue(run.call_args.args[0][3].endswith(".ad.nii.gz"))
        self.assertIn("1 succeeded, 0 skipped, 0 failed", stdout)

    def test_invalid_metric_is_rejected_before_extraction(self):
        self.session()
        for metric in ("unknown", "qa", "fa0", ""):
            with self.subTest(metric=metric), contextlib.redirect_stderr(io.StringIO()) as stderr:
                with patch.object(batch.subprocess, "run") as run:
                    with self.assertRaises(SystemExit) as error:
                        batch.main(["-i", str(self.root), "-r", "annosplit", "-m", metric])
                    self.assertEqual(error.exception.code, 2)
                    run.assert_not_called()
                self.assertIn("invalid choice", stderr.getvalue())

    def test_missing_inputs_are_reported_and_other_sessions_continue(self):
        missing_roi = self.session("a", "1", roi=None)
        missing_dsi = self.session("b", "1", metrics=())
        (missing_dsi / "DSI_studio").rmdir()
        self.session("c", "1", metrics=())
        good = self.session("d", "1")
        status, stdout, stderr, run = self.run_mocked()
        self.assertEqual(status, 0)
        self.assertEqual(run.call_count, 4)
        self.assertIn(str(missing_roi.relative_to(self.root)), stderr)
        self.assertIn("No ROI matching", stderr)
        self.assertIn("DSI_studio folder missing", stderr)
        self.assertIn("[c/1/dwi] No FA map found", stderr)
        self.assertIn("a / 1 (fully skipped: ROI", stdout)
        self.assertIn("b / 1 (fully skipped: DSI_studio folder missing)", stdout)
        self.assertIn("c / 1 (fully skipped: FA, MD, RD, AD missing)", stdout)
        self.assertIn("Summary: 4 succeeded, 6 skipped, 0 failed.", stdout)
        self.assertTrue(all(call.kwargs["cwd"] == good for call in run.call_args_list))

    def test_multiple_roi_matches_fail_session_without_choosing_first(self):
        dwi = self.session("a", "1")
        duplicate = dwi / "other_AnnoSplit_parental.nii"
        duplicate.touch()
        good = self.session("b", "1", metrics=("fa",))
        status, stdout, stderr, run = self.run_mocked()
        self.assertEqual(status, 1)
        self.assertIn("Multiple matching ROIs", stderr)
        self.assertIn(duplicate.name, stderr)
        self.assertIn("scan_AnnoSplit_parental.nii.gz", stderr)
        self.assertIn("  a / 1\n", stdout)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(run.call_args.kwargs["cwd"], good)
        self.assertIn("1 succeeded, 3 skipped, 1 failed", stdout)

    def test_output_collisions_skip_all_conflicting_maps_but_allow_others(self):
        dwi = self.session(metrics=("fa", "md"))
        (dwi / "DSI_studio" / "scan.fz.fa.nii").touch()
        output = dwi / "DSI_studio" / "scan.fz.fa_AnnoSplit_parental.csv"
        output.write_bytes(b"existing CSV")
        status, stdout, stderr, run = self.run_mocked()
        self.assertEqual(status, 1)
        self.assertIn("Colliding CSV outputs", stderr)
        self.assertEqual(run.call_count, 1)
        self.assertTrue(run.call_args.args[0][3].endswith(".md.nii.gz"))
        self.assertEqual(output.read_bytes(), b"existing CSV")
        self.assertIn("1 succeeded, 2 skipped, 2 failed", stdout)

    def test_child_failure_and_launch_error_do_not_stop_remaining_extractions(self):
        self.session("a", "1", metrics=("fa", "md"))
        self.session("b", "1", metrics=("fa",))
        outcomes = [SimpleNamespace(returncode=7), OSError("launch failed"), SimpleNamespace(returncode=0)]
        status, stdout, stderr, run = self.run_mocked(side_effect=outcomes)
        self.assertEqual(status, 1)
        self.assertEqual(run.call_count, 3)
        self.assertIn("exit code 7", stderr)
        self.assertIn("launch failed", stderr)
        self.assertEqual(stdout.count("  a / 1\n"), 1)
        failed_block = stdout.split("Failed sessions:\n")[1].split("Skipped sessions:\n")[0]
        self.assertNotIn("  b / 1", failed_block)
        self.assertIn("a / 1 (partially skipped: RD, AD missing)", stdout)
        self.assertIn("b / 1 (partially skipped: MD, RD, AD missing)", stdout)
        self.assertIn("1 succeeded, 5 skipped, 2 failed", stdout)

    def test_empty_cohort_or_no_success_returns_failure(self):
        status, stdout, stderr, run = self.run_mocked()
        self.assertEqual(status, 1)
        self.assertIn("No subject/session/dwi folders", stdout)
        self.assertIn("0 succeeded", stdout)
        run.assert_not_called()
        self.session(metrics=())
        status, _, _, run = self.run_mocked()
        self.assertEqual(status, 1)
        run.assert_not_called()

    def test_summary_lists_full_and_partial_skips_with_failures_once_and_sorted(self):
        self.session("b", "1", metrics=("fa",))
        self.session("a", "1", roi=None)
        self.session("a", "2", metrics=())
        self.session("c", "1")
        failing = self.session("d", "1", metrics=("fa",))

        def outcome(command, **kwargs):
            return SimpleNamespace(returncode=1 if kwargs["cwd"] == failing else 0)

        status, stdout, stderr, run = self.run_mocked(side_effect=outcome)
        self.assertEqual(status, 1)
        self.assertEqual(run.call_count, 6)
        self.assertIn("No ROI matching", stderr)
        self.assertIn("Extraction failed", stderr)
        self.assertEqual(stdout,
                         "Summary: 5 succeeded, 11 skipped, 1 failed.\n"
                         "Failed sessions:\n  d / 1\n"
                         "Skipped sessions:\n"
                         "  a / 1 (fully skipped: ROI matching '*annosplit_parental.nii[.gz]' missing)\n"
                         "  a / 2 (fully skipped: FA, MD, RD, AD missing)\n"
                         "  b / 1 (partially skipped: MD, RD, AD missing)\n"
                         "  d / 1 (partially skipped: MD, RD, AD missing)\n")

    def test_invalid_inputs_fail_before_running_any_extractor(self):
        self.session()
        for args, message in (
            (["-i", str(self.folder / "absent"), "-r", "annosplit"], "Proc folder"),
            (["-i", str(self.translator), "-r", "annosplit"], "Proc folder"),
            (["-i", str(self.root), "-r", " "], "ROI name must not be empty"),
            (["-i", str(self.root), "-r", "annosplit", "-t", str(self.folder / "absent")], "Required file"),
            (["-i", str(self.root), "-r", "annosplit", "-t", str(self.folder)], "Required file"),
        ):
            with self.subTest(args=args), contextlib.redirect_stderr(io.StringIO()) as stderr, \
                    contextlib.redirect_stdout(io.StringIO()) as stdout:
                with patch.object(batch.subprocess, "run") as run:
                    self.assertEqual(batch.main(args), 1)
                    run.assert_not_called()
                self.assertEqual(stderr.getvalue(), "")
                self.assertIn(message, stdout.getvalue())
                self.assertTrue(stdout.getvalue().startswith("Summary:"))

    def test_required_arguments(self):
        for args in ([], ["-i", str(self.root)], ["-r", "annosplit"]):
            with self.subTest(args=args), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    batch.main(args)
                self.assertEqual(error.exception.code, 2)

    def test_end_to_end_exports_each_session_and_preserves_roi_choice(self):
        git_information = "AIDAmri test Git information\nCommit: cohort-test\n"
        (self.root / "AIDAmri_git_information.txt").write_text(git_information, encoding="utf-8")
        sessions = [self.session("sub-A", "ses-1"), self.session("sub-A", "ses-2"),
                    self.session("sub-B", "ses-1")]
        labels = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        for dwi in sessions:
            self.log_path(dwi).write_text("Existing pipeline log\n")
            nib.save(nib.Nifti1Image(labels, np.eye(4)), dwi / "scan_AnnoSplit_parental.nii.gz")
            # A detailed atlas must not be used for the parental request.
            (dwi / "wrong_AnnoSplit.nii.gz").write_bytes(b"not a NIfTI")
            for index, metric in enumerate(batch.METRICS, start=1):
                values = np.array([np.nan, index, index + 2, index + 4], dtype=np.float32).reshape(2, 2, 1)
                nib.save(nib.Nifti1Image(values, np.eye(4)), dwi / "DSI_studio" / f"scan.fz.{metric}.nii.gz")
        result = subprocess.run(
            [sys.executable, "-B", str(SCRIPT), "--input", self.root.name,
             "--roi_name", "annosplit_parental", "--translatorTXT", self.translator.name],
            cwd=self.folder, capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stderr, "")
        self.assertEqual(result.stdout, "Summary: 12 succeeded, 0 skipped, 0 failed.\n")
        for dwi in sessions:
            for index, metric in enumerate(batch.METRICS, start=1):
                output = dwi / "DSI_studio" / f"scan.fz.{metric}_AnnoSplit_parental.csv"
                with output.open(encoding="utf-8", newline="") as stream:
                    self.assertEqual(list(csv.reader(stream)), [
                        ["roi_id", "acronym", "median"],
                        ["22", "L_PTLp", str(float(index + 1))],
                        ["2022", "R_PTLp", str(float(index + 4))],
                    ])
            self.assertEqual(len(list((dwi / "DSI_studio").glob("*.csv"))), 4)
            log = self.log_text(dwi)
            self.assertTrue(log.startswith("Existing pipeline log\n"))
            self.assertEqual(log.count("=== DTIdata_extract.py "), 1)
            self.assertEqual(log.count("Done:"), 4)
            self.assertIn(str(self.translator), log)
            self.assertNotIn("batch_DTIdata_extract.py", log)
            self.assertEqual(log.count(git_information), 1)
            self.assertTrue(log.endswith(git_information))
            manifest_file = dwi / f".{dwi.parent.parent.name}_{dwi.parent.name}_dwi_aidamri_manifest.json"
            manifest = json.loads(manifest_file.read_text())
            stage = manifest["stages"]["extraction"]
            self.assertEqual(stage["created_files"], [
                f"DSI_studio/scan.fz.{metric}_AnnoSplit_parental.csv" for metric in sorted(batch.METRICS)
            ])
            self.assertIn(self.log_path(dwi).name, stage["modified_files"])

    def test_end_to_end_selected_metric_exports_only_fa(self):
        dwi = self.session("sub-A", "ses-1", metrics=(*batch.METRICS, "fa0"))
        labels = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        values = np.array([np.nan, 1, 3, 5], dtype=np.float32).reshape(2, 2, 1)
        nib.save(nib.Nifti1Image(labels, np.eye(4)), dwi / "scan_AnnoSplit_parental.nii.gz")
        nib.save(nib.Nifti1Image(values, np.eye(4)), dwi / "DSI_studio" / "scan.fz.fa.nii.gz")
        # Other metrics remain invalid placeholders and must never be loaded.
        result = subprocess.run(
            [sys.executable, "-B", str(SCRIPT), "-i", str(self.root),
             "-r", "annosplit_parental", "-t", str(self.translator), "-m", "FA"],
            cwd=self.folder, capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stderr, "")
        self.assertEqual(result.stdout, "Summary: 1 succeeded, 0 skipped, 0 failed.\n")
        self.assertEqual(
            [path.name for path in (dwi / "DSI_studio").glob("*.csv")],
            ["scan.fz.fa_AnnoSplit_parental.csv"],
        )

    def test_one_header_per_batch_run_after_first_metric_fails(self):
        dwi = self.session("sub-A", "ses-1", metrics=("fa", "md"))
        labels = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        values = np.array([np.nan, 1, 3, 5], dtype=np.float32).reshape(2, 2, 1)
        nib.save(nib.Nifti1Image(labels, np.eye(4)), dwi / "scan_AnnoSplit_parental.nii.gz")
        # The first metric remains an invalid placeholder; MD must still succeed.
        nib.save(nib.Nifti1Image(values, np.eye(4)), dwi / "DSI_studio" / "scan.fz.md.nii.gz")
        self.log_path(dwi).write_text("Existing pipeline entries\n")
        for run_number in (1, 2):
            result = subprocess.run(
                [sys.executable, "-B", str(SCRIPT), "-i", str(self.root), "-r", "annosplit_parental"],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 1)
            self.assertIn("1 succeeded, 2 skipped, 1 failed", result.stdout)
            log = self.log_text(dwi)
            self.assertTrue(log.startswith("Existing pipeline entries\n"))
            self.assertEqual(log.count("=== DTIdata_extract.py "), run_number)
            self.assertEqual(log.count("Image:"), 2 * run_number)
            self.assertEqual(log.count("Done:"), run_number)
            self.assertLess(log.index("=== DTIdata_extract.py "), log.index("Image:"))

    def test_git_footer_write_failure_marks_the_session_failed_once(self):
        dwi = self.session()
        with patch.object(batch, "append_git_information", side_effect=OSError("Git append denied")):
            status, stdout, stderr, run = self.run_mocked()
        self.assertEqual(status, 1)
        self.assertEqual(run.call_count, 4)
        self.assertIn("4 succeeded, 0 skipped, 1 failed", stdout)
        self.assertEqual(stderr.count("Cannot append Git information"), 1)
        self.assertEqual(self.log_text(dwi).count("=== DTIdata_extract.py "), 1)

    def test_unwritable_header_fails_maps_and_other_sessions_continue(self):
        bad = self.session("a", "1")
        (bad / "analysis.log").mkdir()
        good = self.session("b", "1")
        status, stdout, stderr, run = self.run_mocked()
        self.assertEqual(status, 1)
        self.assertEqual(stderr.count("Cannot write session header"), 4)
        self.assertIn("[a/1/dwi]", stderr)
        self.assertEqual(run.call_count, 4)
        self.assertEqual({call.kwargs["cwd"] for call in run.call_args_list}, {good})
        self.assertEqual(stdout, "Summary: 4 succeeded, 0 skipped, 4 failed.\nFailed sessions:\n  a / 1\n")
        self.assertTrue((bad / "analysis.log").is_dir())
        self.assertEqual(self.log_text(good).count("=== DTIdata_extract.py "), 1)

    def test_child_errors_before_logging_are_visible_with_session_header(self):
        dwi = self.session(metrics=("fa",))
        extractor = self.folder / "DTIdata_extract.py"
        extractor.write_text(
            "import sys\nprint('early stdout')\nprint('early stderr', file=sys.stderr)\nsys.exit(3)\n",
            encoding="utf-8",
        )
        driver = self.folder / "batch_DTIdata_extract.py"
        driver.write_text(SCRIPT.read_text(encoding="utf-8"), encoding="utf-8")
        result = subprocess.run(
            [sys.executable, "-B", str(driver), "-i", str(self.root), "-r", "annosplit_parental", "-m", "fa"],
            capture_output=True, text=True, check=False, env=self.copied_driver_env(),
        )
        self.assertEqual(result.returncode, 1)
        self.assertEqual(result.stdout, "Summary: 0 succeeded, 0 skipped, 1 failed.\n"
                                        "Failed sessions:\n  mouse A / PT 3\n")
        self.assertEqual(result.stderr.count("early stderr"), 1)
        self.assertNotIn("early stdout", result.stdout)
        self.assertIn("exit code 3", result.stderr)
        self.assertEqual(self.log_text(dwi).count("=== DTIdata_extract.py "), 1)
        self.assertNotIn("early stdout", self.log_text(dwi))

    def test_mixed_real_cohort_logs_errors_and_lists_failed_session_only_once(self):
        git_information = "AIDAmri Git information for mixed cohort\n"
        (self.root / "AIDAmri_git_information.txt").write_text(git_information, encoding="utf-8")
        bad = self.session("sub-bad", "ses-1", metrics=("fa",))
        good = self.session("sub-good", "ses-2", metrics=("fa",))
        # Both invalid maps in the first session must be reported as one failed session.
        (bad / "DSI_studio" / "other.fa.nii").touch()
        labels = np.array([0, 22, 2022, 9999], dtype=np.int16).reshape(2, 2, 1)
        values = np.array([np.nan, 1, 3, 5], dtype=np.float32).reshape(2, 2, 1)
        nib.save(nib.Nifti1Image(labels, np.eye(4)), good / "scan_AnnoSplit_parental.nii.gz")
        nib.save(nib.Nifti1Image(values, np.eye(4)), good / "DSI_studio" / "scan.fz.fa.nii.gz")
        result = subprocess.run(
            [sys.executable, "-B", str(SCRIPT), "-i", str(self.root), "-r", "annosplit_parental",
             "-t", str(self.translator), "-m", "fa"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("Error:", result.stderr)
        self.assertEqual(result.stderr.count("No acronym for ROI ID 9999"), 1)
        self.assertEqual(result.stdout, "Summary: 1 succeeded, 0 skipped, 2 failed.\n"
                                        "Failed sessions:\n  sub-bad / ses-1\n")
        self.assertEqual(self.log_text(bad).count("=== DTIdata_extract.py "), 1)
        self.assertIn("Error:", self.log_text(bad))
        self.assertEqual(self.log_text(good).count("No acronym for ROI ID 9999"), 1)
        self.assertEqual(self.log_text(good).count("Done:"), 1)
        self.assertTrue((good / "sub-good_ses-2_dwi_analysis.log").is_file())
        self.assertNotIn("batch_DTIdata_extract.py", self.log_text(bad))
        self.assertNotIn("batch_DTIdata_extract.py", self.log_text(good))
        for dwi in (bad, good):
            self.assertEqual(self.log_text(dwi).count(git_information), 1)
            self.assertTrue(self.log_text(dwi).endswith(git_information))

    def test_child_environment_is_scoped_to_the_tracked_folder(self):
        untracked = self.session(metrics=("fa",))
        tracked = self.session("sub-A", "ses-1", metrics=("fa",))
        with patch.dict(os.environ, {batch.PARENT_TRACKING_ENV: "inherited-marker",
                                    batch.PARENT_GIT_LOG_ENV: "foreign-session",
                                    batch.PARENT_HEADER_LOG_ENV: "foreign-session", "DTI_TEST_ENV": "keep"}):
            before = os.environ.copy()
            status, _, stderr, run = self.run_mocked(extra=("-m", "fa"))
            self.assertEqual(status, 0, stderr)
            self.assertEqual(dict(os.environ), before)
            self.assertEqual(run.call_count, 2)
            for call in run.call_args_list:
                env = call.kwargs["env"]
                self.assertEqual(env["DTI_TEST_ENV"], "keep")
                self.assertEqual(env[batch.PARENT_GIT_LOG_ENV], str(call.kwargs["cwd"].resolve()))
                self.assertEqual(env[batch.PARENT_HEADER_LOG_ENV], str(call.kwargs["cwd"].resolve()))
                if call.kwargs["cwd"] == tracked:
                    self.assertEqual(env[batch.PARENT_TRACKING_ENV], str(tracked.resolve()))
                else:
                    self.assertEqual(call.kwargs["cwd"], untracked)
                    self.assertNotIn(batch.PARENT_TRACKING_ENV, env)

    def test_fully_skipped_sessions_do_not_receive_git_blocks_even_with_existing_logs(self):
        missing_roi = self.session("a", "1", roi=None)
        missing_dsi = self.session("b", "1", metrics=())
        (missing_dsi / "DSI_studio").rmdir()
        missing_metric = self.session("c", "1", metrics=())
        colliding = self.session("d", "1", metrics=("fa",))
        (colliding / "DSI_studio" / "scan.fz.fa.nii").touch()
        for dwi in (missing_roi, missing_dsi, missing_metric, colliding):
            self.log_path(dwi).write_text("Existing pipeline log\n")
        with patch.object(batch, "append_git_information") as append:
            status, _, _, run = self.run_mocked(extra=("-m", "fa"))
        self.assertEqual(status, 1)
        run.assert_not_called()
        append.assert_not_called()
        for dwi in (missing_roi, missing_dsi, missing_metric, colliding):
            self.assertEqual(self.log_text(dwi), "Existing pipeline log\n")

    def test_git_block_is_written_after_all_attempts_even_on_interrupt(self):
        dwi = self.session("sub-A", "ses-1", metrics=("fa", "md"))
        git_information = "Git information after interrupted session\n"
        source = self.root / "AIDAmri_git_information.txt"
        source.write_text(git_information)
        self.log_path(dwi).write_text("Existing pipeline log\n")
        calls = 0

        def interrupted_child(*args, **kwargs):
            nonlocal calls
            calls += 1
            # No Git block is written between metrics.
            self.assertNotIn(git_information, self.log_text(dwi))
            with self.log_path(dwi).open("a") as log:
                log.write(f"child attempt {calls}\n")
            if calls == 2:
                raise KeyboardInterrupt
            return SimpleNamespace(returncode=0)

        with patch("common.script_logging.BUILD_GIT_INFO_PATH", source), \
                patch.object(batch.subprocess, "run", side_effect=interrupted_child), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(KeyboardInterrupt):
                batch.main(["-i", str(self.root), "-r", "annosplit_parental"])
        self.assertEqual(calls, 2)
        self.assertEqual(self.log_text(dwi).count(git_information), 1)
        self.assertTrue(self.log_text(dwi).endswith(git_information))
        self.assertEqual(self.log_text(dwi).count("=== DTIdata_extract.py "), 1)

    def test_exactly_one_parent_tracker_per_extraction_and_no_child_tracker(self):
        dwi = self.session("sub-A", "ses-1", metrics=("fa", "ad"))
        labels = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        values = np.array([np.nan, 1, 3, 5], dtype=np.float32).reshape(2, 2, 1)
        nib.save(nib.Nifti1Image(labels, np.eye(4)), dwi / "scan_AnnoSplit_parental.nii.gz")
        for metric in ("fa", "ad"):
            nib.save(nib.Nifti1Image(values, np.eye(4)), dwi / "DSI_studio" / f"scan.fz.{metric}.nii.gz")
        trackers = []
        original_start = batch.start_output_tracking

        def start(*args):
            tracker = original_start(*args)
            trackers.append(tracker)
            return tracker

        def inline_child(command, **kwargs):
            self.assertEqual(kwargs["env"][extract.PARENT_TRACKING_ENV], str(dwi))
            with patch.dict(os.environ, kwargs["env"], clear=True), contextlib.redirect_stdout(io.StringIO()):
                return SimpleNamespace(returncode=extract.main(command[2:]))

        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(batch, "start_output_tracking", side_effect=start) as parent_start, \
                patch.object(extract, "start_output_tracking", side_effect=AssertionError("double tracking")) as child_start, \
                patch.object(batch.subprocess, "run", side_effect=inline_child), \
                contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = batch.main(["-i", str(self.root), "-r", "annosplit_parental", "-t", str(self.translator)])
        self.assertEqual(status, 0, stderr.getvalue())
        self.assertEqual(parent_start.call_count, 2)
        child_start.assert_not_called()
        self.assertTrue(all(tracker._finalized for tracker in trackers))
        self.assertIn("2 succeeded, 2 skipped, 0 failed", stdout.getvalue())
        stage = self.manifest(dwi)["stages"]["extraction"]
        self.assertEqual(stage["created_files"], [
            "DSI_studio/scan.fz.ad_AnnoSplit_parental.csv", "DSI_studio/scan.fz.fa_AnnoSplit_parental.csv",
            "sub-A_ses-1_dwi_analysis.log",
        ])
        self.assertEqual(self.log_text(dwi).count("Manifest tracking: parent batch"), 2)
        self.assertEqual(self.log_text(dwi).count("=== DTIdata_extract.py "), 1)

    def test_tracker_start_failure_prevents_child_and_counts_once(self):
        self.session("sub-A", "ses-1", metrics=("fa",))
        with patch.object(batch, "start_output_tracking", side_effect=OSError("tracking denied")) as start:
            status, stdout, stderr, run = self.run_mocked(extra=("-m", "fa"))
        self.assertEqual(status, 1)
        self.assertEqual(start.call_count, 1)
        run.assert_not_called()
        self.assertIn("Cannot start manifest tracking", stderr)
        self.assertIn("tracking denied", stderr)
        self.assertEqual(stdout, "Summary: 0 succeeded, 0 skipped, 1 failed.\nFailed sessions:\n  sub-A / ses-1\n")

    def test_finalization_failure_counts_once_even_if_child_also_fails(self):
        self.session("sub-A", "ses-1", metrics=("fa",))
        for outcome in (SimpleNamespace(returncode=0), SimpleNamespace(returncode=3), OSError("launch denied")):
            tracker = SimpleNamespace(finalize=Mock(side_effect=OSError("manifest denied")))
            with self.subTest(outcome=outcome), patch.object(batch, "start_output_tracking", return_value=tracker):
                status, stdout, stderr, run = self.run_mocked(extra=("-m", "fa"), side_effect=[outcome])
                self.assertEqual(status, 1)
                self.assertEqual(run.call_count, 1)
                tracker.finalize.assert_called_once_with()
                self.assertIn("Cannot finalize extraction manifest", stderr)
                self.assertIn("manifest denied", stderr)
                self.assertIn("Summary: 0 succeeded, 0 skipped, 1 failed.", stdout)

    def test_tracker_finalizes_on_keyboard_interrupt(self):
        self.session("sub-A", "ses-1", metrics=("fa",))
        tracker = SimpleNamespace(finalize=Mock())
        with patch.object(batch, "start_output_tracking", return_value=tracker), \
                patch.object(batch.subprocess, "run", side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):
                batch.main(["-i", str(self.root), "-r", "annosplit_parental", "-m", "fa"])
        tracker.finalize.assert_called_once_with()

    @unittest.skipUnless(hasattr(signal, "SIGKILL"), "SIGKILL requires a Unix platform")
    def test_killed_child_outputs_are_tracked_by_surviving_batch(self):
        dwi = self.session("sub-killed", "ses-1", metrics=("fa",))
        driver = self.folder / "batch_DTIdata_extract.py"
        driver.write_text(SCRIPT.read_text(encoding="utf-8"), encoding="utf-8")
        (self.folder / "DTIdata_extract.py").write_text(
            "import os, signal\nfrom pathlib import Path\n"
            "Path('DSI_studio/partial.csv').write_text('incomplete output')\n"
            "os.kill(os.getpid(), signal.SIGKILL)\n", encoding="utf-8",
        )
        result = subprocess.run(
            [sys.executable, "-B", str(driver), "-i", str(self.root), "-r", "annosplit_parental", "-m", "fa"],
            capture_output=True, text=True, check=False, env=self.copied_driver_env(),
        )
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(f"exit code {-signal.SIGKILL}", result.stderr)
        self.assertEqual(result.stdout, "Summary: 0 succeeded, 0 skipped, 1 failed.\n"
                                        "Failed sessions:\n  sub-killed / ses-1\n")
        self.assertEqual(self.manifest(dwi)["stages"]["extraction"]["created_files"], [
            "DSI_studio/partial.csv", self.log_path(dwi).name,
        ])
        self.assertEqual(self.log_text(dwi).count("=== DTIdata_extract.py "), 1)


if __name__ == "__main__":
    unittest.main()
