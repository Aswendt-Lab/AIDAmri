"""Regression tests for AIDAmri manifests and the manifest-based reset.

This module is used for development and quality assurance; it is not required
to process MRI data. It verifies that newly created and modified files are
recorded in the modality-specific manifest and that ``Reset_proc_folder.py``
deletes only registered outputs from stages after the selected target phase.

The tests also cover manifest naming, T2map support, safe restoration of base
NIfTI files, preservation of unmanaged files, and missing manifests inside or
outside the selected modality. Every test uses an automatically created
temporary directory; real project data is never read or modified.

Run from the repository root::

    python -m unittest discover -s tests -v
"""

import contextlib
import io
import stat
import sys
import tempfile
import unittest
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from common.artifact_manifest import (
    LOCK_FILENAME,
    OutputTracker,
    load_manifest,
    manifest_filename,
)
from helper_tools.Reset_proc_folder import (
    apply_reset_plan,
    build_reset_plan,
    main as reset_main,
    print_issue_summary,
    print_reset_plan,
    validate_project_manifests,
)


class ArtifactManifestTests(unittest.TestCase):
    def test_tracker_records_created_and_modified_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            folder = Path(temp_dir) / "sub-test" / "ses-test" / "anat"
            folder.mkdir(parents=True)
            existing = folder / "input.nii.gz"
            existing.write_text("before", encoding="utf-8")

            tracker = OutputTracker.start(folder, "anat", "preprocessing")
            existing.write_text("after", encoding="utf-8")
            (folder / "output.nii.gz").write_text("new", encoding="utf-8")
            tracker.finalize()

            self.assertEqual(
                manifest_filename(folder, "anat"),
                ".sub-test_ses-test_anat_aidamri_manifest.json",
            )
            manifest_file = folder / manifest_filename(folder, "anat")
            self.assertEqual(
                stat.S_IMODE(manifest_file.stat().st_mode),
                0o644,
            )
            stage = load_manifest(folder, "anat")["stages"]["preprocessing"]
            self.assertEqual(stage["created_files"], ["output.nii.gz"])
            self.assertEqual(stage["modified_files"], ["input.nii.gz"])

    def test_reset_deletes_only_manifested_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "anat"
            (folder / "brkraw").mkdir(parents=True)
            raw = folder / "input_T2w.nii.gz"
            raw.write_text("raw", encoding="utf-8")

            preprocessing = OutputTracker.start(folder, "anat", "preprocessing")
            preprocessing_output = folder / "preprocessed.nii.gz"
            preprocessing_output.write_text("pre", encoding="utf-8")
            preprocessing.finalize()

            registration = OutputTracker.start(folder, "anat", "registration")
            registration_dir = folder / "registration_outputs"
            registration_dir.mkdir()
            registration_output = registration_dir / "registered.nii.gz"
            registration_output.write_text("reg", encoding="utf-8")
            raw.write_text("changed by registration", encoding="utf-8")
            registration.finalize()

            unknown = registration_dir / "notes.txt"
            unknown.write_text("user data", encoding="utf-8")
            muted_sidecars = [
                folder / "sub-01_ses-01_T2w.json",
                folder / "sub-01_ses-01_EPI.json",
                folder / "sub-01_ses-01_dwi.json",
                folder / "sub-01_ses-01_dwi.bvec",
                folder / "sub-01_ses-01_dwi.bval",
            ]
            for sidecar in muted_sidecars:
                sidecar.write_text("{}", encoding="utf-8")
            stroke_mask = folder / "sub-01_ses-01_BetStroke_mask.nii.gz"
            stroke_mask.write_text("stroke", encoding="utf-8")

            plans = build_reset_plan(project, "anat", "preprocessing")
            self.assertIn(unknown, plans[0]["unknown_files"])
            self.assertIn(stroke_mask, plans[0]["unknown_files"])
            for sidecar in muted_sidecars:
                self.assertIn(sidecar, plans[0]["unknown_files"])
            report = io.StringIO()
            with contextlib.redirect_stdout(report):
                print_reset_plan(plans, "anat", "preprocessing")
            report_text = report.getvalue()
            self.assertIn(
                f"Not managed by the manifest; preserved: {unknown}",
                report_text,
            )
            self.assertIn(
                "[INFO] Stroke mask files are present in 1 selected folder.",
                report_text,
            )
            self.assertIn(
                f"Not managed by the manifest; preserved: {stroke_mask}",
                report_text,
            )
            for sidecar in muted_sidecars:
                self.assertIn(
                    f"Not managed by the manifest; preserved: {sidecar}",
                    report_text,
                )
            summary = io.StringIO()
            with contextlib.redirect_stdout(summary):
                print_issue_summary(plans)
            summary_text = summary.getvalue()
            self.assertIn(
                "[INFO] Stroke mask files are present in 1 selected folder "
                "and were preserved.",
                summary_text,
            )
            self.assertIn("sub-01 | ses-01 | anat", summary_text)
            self.assertIn(
                f"[INFO] Not managed by the manifest; review manually: {unknown}",
                summary_text,
            )
            self.assertIn(
                f"[WARNING] Modified by registration; original content cannot be "
                f"restored: {raw}",
                summary_text,
            )
            for sidecar in muted_sidecars:
                self.assertNotIn(str(sidecar), summary_text)
            self.assertNotIn(str(stroke_mask), summary_text)
            with contextlib.redirect_stdout(io.StringIO()):
                apply_reset_plan(plans)

            self.assertTrue(raw.exists())
            self.assertTrue(preprocessing_output.exists())
            self.assertFalse(registration_output.exists())
            self.assertTrue(unknown.exists())
            for sidecar in muted_sidecars:
                self.assertTrue(sidecar.exists())
            self.assertTrue(stroke_mask.exists())
            self.assertTrue(registration_dir.exists())
            self.assertIsNotNone(load_manifest(folder, "anat"))
            self.assertTrue((folder / LOCK_FILENAME).exists())

    def test_folder_without_manifest_is_skipped(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "anat"
            (folder / "brkraw").mkdir(parents=True)
            unknown = folder / "unknown.nii.gz"
            unknown.write_text("keep", encoding="utf-8")

            plans = build_reset_plan(project, "anat", "base")

            self.assertTrue(plans[0]["skip"])
            self.assertTrue(unknown.exists())

    def test_delete_unmanaged_requires_opt_in_and_removes_unknown_paths(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "anat"
            folder.mkdir(parents=True)
            raw = folder / "input_T2w.nii.gz"
            raw.write_text("original", encoding="utf-8")

            tracker = OutputTracker.start(folder, "anat", "preprocessing")
            brkraw = folder / "brkraw"
            brkraw.mkdir()
            backup = brkraw / raw.name
            backup.write_text("original", encoding="utf-8")
            raw.write_text("preprocessed", encoding="utf-8")
            known_output = folder / "known_output.nii.gz"
            known_output.write_text("known", encoding="utf-8")
            tracker.finalize()

            sidecar = folder / "sub-01_ses-01_T2w.json"
            sidecar.write_text("{}", encoding="utf-8")
            bval = folder / "sub-01_ses-01_dwi.bval"
            bval.write_text("0", encoding="utf-8")
            stroke_mask = folder / "sub-01_ses-01_BetStroke_mask.nii.gz"
            stroke_mask.write_text("stroke", encoding="utf-8")
            manual_mask = folder / "manual_mask.nii.gz"
            manual_mask.write_text("manual", encoding="utf-8")
            nested_directory = folder / "custom" / "nested"
            nested_directory.mkdir(parents=True)
            note = nested_directory / "notes.txt"
            note.write_text("manual", encoding="utf-8")
            protected_directory = folder / "protected" / "nested"
            protected_directory.mkdir(parents=True)
            bvec = protected_directory / "sub-01_ses-01_dwi.bvec"
            bvec.write_text("1 0 0", encoding="utf-8")
            empty_directory = folder / "empty_custom_directory"
            empty_directory.mkdir()

            safe_plan = build_reset_plan(project, "anat", "processing")[0]
            self.assertFalse(safe_plan["delete_unmanaged"])
            self.assertIn(sidecar, safe_plan["unknown_files"])
            self.assertNotIn(sidecar, safe_plan["files"])

            dry_run_output = io.StringIO()
            with contextlib.redirect_stdout(dry_run_output):
                exit_code = reset_main(
                    [
                        "--input",
                        str(project),
                        "--mode",
                        "anat",
                        "--phase",
                        "processing",
                        "--delete-unmanaged",
                        "--dry-run",
                    ]
                )
            self.assertEqual(exit_code, 0)
            dry_run_text = dry_run_output.getvalue()
            self.assertIn("[DANGER] --delete-unmanaged is active", dry_run_text)
            self.assertIn(
                f"[INFO] Protected from --delete-unmanaged; preserved: {sidecar}",
                dry_run_text,
            )
            self.assertIn(
                f"[INFO] Protected from --delete-unmanaged; preserved: {stroke_mask}",
                dry_run_text,
            )
            self.assertTrue(sidecar.exists())

            plans = build_reset_plan(
                project,
                "anat",
                "processing",
                delete_unmanaged=True,
            )
            plan = plans[0]
            self.assertEqual(
                set(plan["unknown_files"]),
                {sidecar, bval, stroke_mask, manual_mask, note, bvec},
            )
            self.assertEqual(
                set(plan["protected_unmanaged_files"]),
                {sidecar, bval, stroke_mask, bvec},
            )
            self.assertEqual(
                set(plan["unmanaged_files_to_delete"]),
                {manual_mask, note},
            )
            self.assertEqual(
                set(plan["unknown_directories"]),
                {
                    folder / "custom",
                    nested_directory,
                    folder / "protected",
                    protected_directory,
                    empty_directory,
                },
            )
            self.assertEqual(
                set(plan["unmanaged_directories_to_delete"]),
                {folder / "custom", nested_directory, empty_directory},
            )

            with contextlib.redirect_stdout(io.StringIO()):
                succeeded = apply_reset_plan(plans)

            self.assertTrue(succeeded)
            self.assertTrue(sidecar.exists())
            self.assertTrue(bval.exists())
            self.assertTrue(stroke_mask.exists())
            self.assertTrue(bvec.exists())
            self.assertTrue(protected_directory.exists())
            self.assertTrue((folder / "protected").exists())
            self.assertFalse(manual_mask.exists())
            self.assertFalse(note.exists())
            self.assertFalse(nested_directory.exists())
            self.assertFalse((folder / "custom").exists())
            self.assertFalse(empty_directory.exists())
            self.assertTrue(raw.exists())
            self.assertTrue(backup.exists())
            self.assertTrue(known_output.exists())
            self.assertIsNotNone(load_manifest(folder, "anat"))

    def test_base_reset_restores_original_nifti_before_deleting_backup(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "anat"
            folder.mkdir(parents=True)
            raw = folder / "input_T2w.nii.gz"
            raw.write_text("original", encoding="utf-8")

            preprocessing = OutputTracker.start(folder, "anat", "preprocessing")
            brkraw = folder / "brkraw"
            brkraw.mkdir()
            backup = brkraw / raw.name
            backup.write_text("original", encoding="utf-8")
            raw.write_text("preprocessed", encoding="utf-8")
            preprocessing_output = folder / "preprocessed.nii.gz"
            preprocessing_output.write_text("pre", encoding="utf-8")
            preprocessing_log = folder / "preprocess.log"
            preprocessing_log.write_text("first run", encoding="utf-8")
            preprocessing.finalize()

            # On later runs a stage-owned log is both created and modified in
            # the accumulated manifest. It is still safe to delete on reset.
            preprocessing = OutputTracker.start(folder, "anat", "preprocessing")
            raw.write_text("preprocessed again", encoding="utf-8")
            preprocessing_log.write_text("second run", encoding="utf-8")
            preprocessing.finalize()

            registration = OutputTracker.start(folder, "anat", "registration")
            registration_output = folder / "registered.nii.gz"
            registration_output.write_text("reg", encoding="utf-8")
            registration.finalize()

            plans = build_reset_plan(project, "anat", "base")
            self.assertEqual(plans[0]["restore_files"], [(backup, raw)])
            self.assertEqual(plans[0]["restore_errors"], [])
            self.assertEqual(plans[0]["modified_files"], [])
            manifest_file = folder / manifest_filename(folder, "anat")
            lock_file = folder / LOCK_FILENAME
            self.assertEqual(
                plans[0]["metadata_files"],
                [manifest_file, lock_file],
            )
            self.assertTrue(manifest_file.exists())
            self.assertTrue(lock_file.exists())

            with contextlib.redirect_stdout(io.StringIO()):
                succeeded = apply_reset_plan(plans)

            self.assertTrue(succeeded)
            self.assertEqual(raw.read_text(encoding="utf-8"), "original")
            self.assertFalse(backup.exists())
            self.assertFalse(brkraw.exists())
            self.assertFalse(preprocessing_output.exists())
            self.assertFalse(preprocessing_log.exists())
            self.assertFalse(registration_output.exists())
            self.assertIsNone(load_manifest(folder, "anat"))
            self.assertFalse(manifest_file.exists())
            self.assertFalse(lock_file.exists())

            # The next processing run starts a fresh manifest and lock.
            tracker = OutputTracker.start(folder, "anat", "preprocessing")
            (folder / "new_output.nii.gz").write_text("new", encoding="utf-8")
            tracker.finalize()
            self.assertIsNotNone(load_manifest(folder, "anat"))
            self.assertTrue(lock_file.exists())

    def test_base_reset_aborts_before_changes_when_backup_is_missing(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "anat"
            folder.mkdir(parents=True)
            raw = folder / "input_T2w.nii.gz"
            raw.write_text("original", encoding="utf-8")

            tracker = OutputTracker.start(folder, "anat", "preprocessing")
            brkraw = folder / "brkraw"
            brkraw.mkdir()
            backup = brkraw / raw.name
            backup.write_text("original", encoding="utf-8")
            raw.write_text("preprocessed", encoding="utf-8")
            output = folder / "preprocessed.nii.gz"
            output.write_text("keep until safe reset", encoding="utf-8")
            tracker.finalize()
            backup.unlink()

            report = io.StringIO()
            with contextlib.redirect_stdout(report):
                exit_code = reset_main(
                    [
                        "--input",
                        str(project),
                        "--mode",
                        "anat",
                        "--phase",
                        "base",
                        "--yes",
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn("Base reset aborted", report.getvalue())
            self.assertEqual(raw.read_text(encoding="utf-8"), "preprocessed")
            self.assertTrue(output.exists())
            self.assertIsNotNone(load_manifest(folder, "anat"))
            self.assertTrue((folder / LOCK_FILENAME).exists())

    def test_t2map_outputs_can_be_tracked(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            folder = Path(temp_dir) / "sub-test" / "ses-test" / "t2map"
            folder.mkdir(parents=True)
            tracker = OutputTracker.start(folder, "t2map", "processing")
            (folder / "t2values.csv").write_text("value", encoding="utf-8")
            tracker.finalize()

            manifest = load_manifest(folder, "t2map")
            self.assertEqual(manifest["mode"], "t2map")
            self.assertEqual(
                manifest["stages"]["processing"]["created_files"],
                ["t2values.csv"],
            )

    def test_reset_ignores_missing_manifest_in_other_modality(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            anat = project / "sub-01" / "ses-01" / "anat"
            dwi = project / "sub-01" / "ses-01" / "dwi"
            (anat / "brkraw").mkdir(parents=True)
            (dwi / "brkraw").mkdir(parents=True)

            tracker = OutputTracker.start(anat, "anat", "preprocessing")
            tracker.finalize()

            tracker = OutputTracker.start(anat, "anat", "registration")
            anat_output = anat / "preprocessed.nii.gz"
            anat_output.write_text("pre", encoding="utf-8")
            tracker.finalize()

            self.assertEqual(validate_project_manifests(project, "anat"), [])
            dwi_errors = validate_project_manifests(project, "dwi")
            self.assertEqual(len(dwi_errors), 1)
            self.assertIn(
                "/dwi/.sub-01_ses-01_dwi_aidamri_manifest.json",
                dwi_errors[0],
            )

            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = reset_main(
                    [
                        "--input",
                        str(project),
                        "--mode",
                        "anat",
                        "--phase",
                        "preprocessing",
                        "--yes",
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertFalse(anat_output.exists())

    def test_reset_aborts_when_selected_modality_manifest_is_missing(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            dwi = project / "sub-01" / "ses-01" / "dwi"
            (dwi / "brkraw").mkdir(parents=True)
            unknown = dwi / "unknown.nii.gz"
            unknown.write_text("keep", encoding="utf-8")

            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = reset_main(
                    [
                        "--input",
                        str(project),
                        "--mode",
                        "dwi",
                        "--phase",
                        "base",
                        "--yes",
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertTrue(unknown.exists())


if __name__ == "__main__":
    unittest.main()
