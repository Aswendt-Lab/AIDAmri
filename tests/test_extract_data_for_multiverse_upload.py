"""Tests for collecting Multiverse results for upload."""

import importlib.util
import tempfile
import unittest
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
MODULE_PATH = BIN_DIR / "Extract_Data_for_multiverse_upload.py"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "extract_data_for_multiverse_upload_test",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ExtractMultiverseUploadTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def _create_file(self, root, relative_path, content=b"test data"):
        path = Path(root) / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        return path

    def test_copies_selected_files_directly_into_bids_modality_folder(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            multiverse_root = Path(temp_dir) / "Multiverse"
            project_root = multiverse_root / "proc_data_good"
            output_dir = multiverse_root / "Extracted_data_for_upload"
            multiverse_output = (
                Path("sub-01/ses-1/func/rs-fMRI_niiData/Multiverse_Output")
            )
            corrected = self._create_file(
                project_root,
                multiverse_output
                / "sub-01_registered_on_SIGMA_template_corrected.nii.gz",
                b"4d",
            )
            corrected_mean = self._create_file(
                project_root,
                multiverse_output
                / (
                    "sub-01_registered_on_SIGMA_template_"
                    "temporal_mean_corrected.nii.gz"
                ),
                b"mean",
            )
            self._create_file(
                project_root,
                "Report/Multiverse/Registration/ignored.nii.gz",
            )

            copied, pattern_counts = self.module.copy_upload_files(project_root)

            self.assertEqual(len(copied), 2)
            self.assertEqual(sum(pattern_counts.values()), 2)
            for source in (corrected, corrected_mean):
                destination = (
                    output_dir / "sub-01/ses-1/func" / source.name
                )
                self.assertTrue(destination.is_file())
                self.assertEqual(destination.read_bytes(), source.read_bytes())
            self.assertFalse(
                (output_dir / "sub-01/ses-1/func/rs-fMRI_niiData").exists()
            )
            self.assertFalse((output_dir / "Report").exists())

    def test_flattens_anatomical_file_to_anat_folder(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project_root = Path(temp_dir) / "proc_data_good"
            source = self._create_file(
                project_root,
                "sub-01/ses-1/anat/nested/results/anatomy_upload.nii.gz",
            )

            relative_path = self.module.upload_relative_path(
                source,
                project_root,
            )

            self.assertEqual(
                relative_path,
                Path("sub-01/ses-1/anat/anatomy_upload.nii.gz"),
            )

    def test_copies_complete_directory_into_modality_folder(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            multiverse_root = Path(temp_dir) / "Multiverse"
            project_root = multiverse_root / "proc_data_good"
            output_root = multiverse_root / "Extracted_data_for_upload"
            source_file = self._create_file(
                project_root,
                "sub-01/ses-1/func/deep/results/regr/nested/matrix.mat",
                b"matrix data",
            )
            self._create_file(
                project_root,
                "Report/regr/ignored.mat",
            )
            self.module.DIRECTORY_PATTERNS = (
                "sub-*/ses-*/func/**/regr",
            )

            copied, pattern_counts = self.module.copy_upload_files(project_root)

            destination_directory = output_root / "sub-01/ses-1/func/regr"
            destination_file = destination_directory / "nested/matrix.mat"
            self.assertIn(destination_directory, copied)
            self.assertEqual(
                destination_file.read_bytes(),
                source_file.read_bytes(),
            )
            self.assertEqual(
                pattern_counts[
                    ("directory", "sub-*/ses-*/func/**/regr")
                ],
                1,
            )
            self.assertFalse((output_root / "Report").exists())

    def test_rejects_two_sources_with_same_flattened_destination(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project_root = Path(temp_dir) / "proc_data_good"
            self._create_file(
                project_root,
                "sub-01/ses-1/func/first/duplicate_upload.nii.gz",
            )
            self._create_file(
                project_root,
                "sub-01/ses-1/func/second/duplicate_upload.nii.gz",
            )
            self.module.FILE_PATTERNS = (
                "sub-*/ses-*/func/**/*_upload.nii.gz",
            )

            with self.assertRaisesRegex(
                FileExistsError,
                "same destination",
            ):
                self.module.copy_upload_files(project_root)

    def test_dry_run_does_not_create_output_directory(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            multiverse_root = Path(temp_dir) / "Multiverse"
            project_root = multiverse_root / "proc_data_good"
            self._create_file(
                project_root,
                (
                    "sub-01/ses-1/func/rs-fMRI_niiData/Multiverse_Output/"
                    "sub-01_registered_on_SIGMA_template_corrected.nii.gz"
                ),
            )

            copied, _ = self.module.copy_upload_files(
                project_root,
                dry_run=True,
            )

            self.assertEqual(len(copied), 1)
            self.assertFalse(
                (multiverse_root / "Extracted_data_for_upload").exists()
            )


if __name__ == "__main__":
    unittest.main()
