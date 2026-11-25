"""Regression tests for AIDAmri manifests and the manifest-based reset.

This module is used for development and quality assurance; it is not required
to process MRI data. It verifies that newly created and modified files are
recorded in the modality-specific manifest and that ``Reset_proc_folder.py``
deletes only registered outputs from stages after the selected target phase.

The tests also cover manifest naming, T2map support, preservation of unmanaged
files, and missing manifests inside or outside the selected modality. Every
test uses an automatically created temporary directory; real project data is
never read or modified.

Run from the repository root::

    python -m unittest discover -s tests -v
"""

import contextlib
import importlib.util
import io
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from common.artifact_manifest import OutputTracker, load_manifest, manifest_filename
from helper_tools.Reset_proc_folder import (
    apply_reset_plan,
    build_reset_plan,
    main as reset_main,
    print_issue_summary,
    print_reset_plan,
    normalize_phase,
    parse_args,
    validate_project_manifests,
)


class ArtifactManifestTests(unittest.TestCase):
    def test_extraction_alias_preserves_existing_processing_entries_and_manifest_version(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            folder = Path(temp_dir) / "sub-01" / "ses-01" / "dwi"
            folder.mkdir(parents=True)
            processing = OutputTracker.start(folder, "dwi", "processing")
            (folder / "legacy.csv").write_text("legacy export")
            processing.finalize()
            before = load_manifest(folder, "dwi")

            extraction = OutputTracker.start(folder, "dwi", " ExTrAcT ")
            (folder / "new.csv").write_text("new export")
            extraction.finalize()
            after = load_manifest(folder, "dwi")
            self.assertEqual(after["version"], before["version"])
            self.assertEqual(after["version"], 1)
            self.assertEqual(after["tracking_started_before_stage"], "processing")
            self.assertEqual(after["stages"]["processing"], before["stages"]["processing"])
            self.assertEqual(after["stages"]["extraction"]["created_files"], ["new.csv"])
            self.assertNotIn("extract", after["stages"])
            self.assertEqual(normalize_phase("extract"), "extraction")

    def test_reset_to_processing_deletes_new_extraction_outputs_and_preserves_modified_legacy_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "dwi"
            (folder / "brkraw").mkdir(parents=True)
            processing = OutputTracker.start(folder, "dwi", "processing")
            image = folder / "processed.nii.gz"
            image.write_text("processing output")
            legacy = folder / "legacy.csv"
            legacy.write_text("legacy export")
            log = folder / "process.log"
            log.write_text("pipeline log")
            processing.finalize()
            previous_stage = load_manifest(folder, "dwi")["stages"]["processing"]

            extraction = OutputTracker.start(folder, "dwi", "extraction")
            output = folder / "new.csv"
            output.write_text("new export")
            legacy.write_text("updated legacy export")
            log.write_text("pipeline log plus extraction")
            extraction.finalize()
            plans = build_reset_plan(project, "dwi", "processing")
            self.assertEqual(plans[0]["selected_stages"], ["extraction"])
            with contextlib.redirect_stdout(io.StringIO()):
                apply_reset_plan(plans)

            self.assertFalse(output.exists())
            self.assertEqual(image.read_text(), "processing output")
            self.assertEqual(legacy.read_text(), "updated legacy export")
            self.assertEqual(log.read_text(), "pipeline log plus extraction")
            manifest = load_manifest(folder, "dwi")
            self.assertEqual(manifest["stages"]["processing"], previous_stage)
            self.assertEqual(manifest["stages"]["extraction"]["created_files"], [])
            self.assertEqual(manifest["stages"]["extraction"]["modified_files"], ["legacy.csv", "process.log"])

    def test_reset_to_extraction_is_accepted_by_cli_and_preserves_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "dwi"
            (folder / "brkraw").mkdir(parents=True)
            tracker = OutputTracker.start(folder, "dwi", "extraction")
            output = folder / "values.csv"
            output.write_text("extraction output")
            tracker.finalize()
            args = ["-i", str(project), "-m", "dwi", "-p", "extraction", "--dry-run"]
            self.assertEqual(parse_args(args).phase, "extraction")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(reset_main(args), 0)
                plans = build_reset_plan(project, "dwi", "extract")
                self.assertEqual(plans[0]["selected_stages"], [])
                apply_reset_plan(plans)
            self.assertEqual(output.read_text(), "extraction output")
            self.assertIn("extraction", load_manifest(folder, "dwi")["stages"])

    def test_reset_to_registration_removes_processing_and_extraction_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "dwi"
            (folder / "brkraw").mkdir(parents=True)
            paths = {}
            for stage in ("registration", "processing", "extraction"):
                tracker = OutputTracker.start(folder, "dwi", stage)
                paths[stage] = folder / f"{stage}.txt"
                paths[stage].write_text(stage)
                tracker.finalize()
            plans = build_reset_plan(project, "dwi", "registration")
            self.assertEqual(set(plans[0]["selected_stages"]), {"processing", "extraction"})
            with contextlib.redirect_stdout(io.StringIO()):
                apply_reset_plan(plans)
            self.assertTrue(paths["registration"].exists())
            self.assertFalse(paths["processing"].exists())
            self.assertFalse(paths["extraction"].exists())

    def test_tracking_starting_at_extraction_retains_late_tracking_warning(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "dwi"
            (folder / "brkraw").mkdir(parents=True)
            existing = folder / "untracked.nii.gz"
            existing.write_text("previous processing output")
            tracker = OutputTracker.start(folder, "dwi", "extraction")
            (folder / "values.csv").write_text("export")
            tracker.finalize()
            self.assertEqual(load_manifest(folder, "dwi")["tracking_started_before_stage"], "extraction")
            plans = build_reset_plan(project, "dwi", "registration")
            self.assertTrue(any("tracking started only before stage 'extraction'" in warning
                                for warning in plans[0]["warnings"]))
            with contextlib.redirect_stdout(io.StringIO()):
                apply_reset_plan(plans)
            self.assertEqual(existing.read_text(), "previous processing output")

    def test_multiverse_outputs_are_tracked_in_separate_stage(self):
        module_path = BIN_DIR / "Create_multiverse_output.py"
        spec = importlib.util.spec_from_file_location(
            "create_multiverse_output_test",
            module_path,
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory() as temp_dir:
            subject = Path(temp_dir) / "sub-01" / "ses-1"
            func_root = subject / "func"
            func_output = func_root / "rs-fMRI_niiData"
            anat = subject / "anat"
            func_output.mkdir(parents=True)
            anat.mkdir()
            epi = func_output / "sub-01_task-rest_bold_EPI_mcf_f.nii.gz"
            epi.write_text("input", encoding="utf-8")

            registered = func_output / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template.nii.gz"
            )
            temporal_mean = func_output / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_temporal_mean.nii.gz"
            )

            def fake_transform(**_kwargs):
                registered.write_text("registered", encoding="utf-8")
                return [str(registered)]

            def fake_temporal_mean(_input_path, output_path):
                Path(output_path).write_text("mean", encoding="utf-8")

            with (
                mock.patch.object(
                    module,
                    "apply_inverse_composite_transformation",
                    side_effect=fake_transform,
                ),
                mock.patch.object(
                    module,
                    "compute_temporal_mean",
                    side_effect=fake_temporal_mean,
                ),
            ):
                module.process_subject(str(subject), "template.nii.gz")

            manifest = load_manifest(func_root, "func")
            self.assertEqual(manifest["mode"], "func")
            self.assertEqual(
                manifest["stages"]["multiverse_output"]["created_files"],
                [
                    registered.relative_to(func_root).as_posix(),
                    temporal_mean.relative_to(func_root).as_posix(),
                ],
            )

    def test_reset_to_processing_deletes_only_multiverse_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "func"
            (folder / "brkraw").mkdir(parents=True)

            processing = OutputTracker.start(folder, "func", "processing")
            processing_output = folder / "processed.nii.gz"
            processing_output.write_text("processed", encoding="utf-8")
            processing.finalize()

            multiverse = OutputTracker.start(
                folder,
                "func",
                "multiverse_output",
            )
            multiverse_output = folder / "multiverse.nii.gz"
            multiverse_output.write_text("multiverse", encoding="utf-8")
            multiverse.finalize()

            plans = build_reset_plan(project, "func", "processing")
            self.assertEqual(plans[0]["selected_stages"], ["multiverse_output"])

            with contextlib.redirect_stdout(io.StringIO()):
                apply_reset_plan(plans)

            self.assertTrue(processing_output.exists())
            self.assertFalse(multiverse_output.exists())
            manifest = load_manifest(folder, "func")
            self.assertIn("processing", manifest["stages"])
            self.assertNotIn("multiverse_output", manifest["stages"])

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

            plans = build_reset_plan(project, "anat", "preprocessing")
            self.assertIn(unknown, plans[0]["unknown_files"])
            report = io.StringIO()
            with contextlib.redirect_stdout(report):
                print_reset_plan(plans, "anat", "preprocessing")
            self.assertIn(
                f"Not managed by the manifest; preserved: {unknown}",
                report.getvalue(),
            )
            summary = io.StringIO()
            with contextlib.redirect_stdout(summary):
                print_issue_summary(plans)
            summary_text = summary.getvalue()
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
            with contextlib.redirect_stdout(io.StringIO()):
                apply_reset_plan(plans)

            self.assertTrue(raw.exists())
            self.assertTrue(preprocessing_output.exists())
            self.assertFalse(registration_output.exists())
            self.assertTrue(unknown.exists())
            self.assertTrue(registration_dir.exists())

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
                        "base",
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
