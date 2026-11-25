"""Tests for corrected Multiverse registration and CC QC reports."""

import sys
import tempfile
import unittest
from pathlib import Path

import nibabel as nib
import numpy as np


HELPER_DIR = Path(__file__).resolve().parents[1] / "bin" / "helper_tools"
sys.path.insert(0, str(HELPER_DIR))

import batch_multiverse_qc_reports as reports


class MultiverseQcReportTests(unittest.TestCase):
    def _create_inputs(self, root, atlas_affine=None):
        root = Path(root)
        output_dir = (
            root
            / "sub-01"
            / "ses-1"
            / "func"
            / "rs-fMRI_niiData"
            / "Multiverse_Output"
        )
        output_dir.mkdir(parents=True)
        mean_path = output_dir / (
            "sub-01_task-rest_bold_EPI_mcf_st_f_"
            "registered_on_SIGMA_template_temporal_mean_corrected.nii.gz"
        )
        atlas_path = root / "SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz"

        affine = np.diag([-0.15, -0.15, 0.15, 1.0])
        mean_data = np.arange(8 * 9 * 10, dtype=np.float32).reshape(8, 9, 10)
        atlas_data = np.zeros((8, 9, 10), dtype=np.int16)
        atlas_data[1:7, 1:8, 1:9] = 1
        atlas_data[3:5, 3:6, 3:7] = 891
        atlas_data[4:6, 4:7, 4:8] = 892

        nib.save(nib.Nifti1Image(mean_data, affine), mean_path)
        nib.save(
            nib.Nifti1Image(
                atlas_data,
                affine if atlas_affine is None else atlas_affine,
            ),
            atlas_path,
        )
        return mean_path, atlas_path

    def test_builds_registration_and_cc_reports(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            mean_path, atlas_path = self._create_inputs(temp_dir)

            registration_html, registration_count = (
                reports.build_multiverse_registration_qc_report(
                    temp_dir,
                    atlas_path=atlas_path,
                    n_slices=2,
                )
            )
            cc_html, cc_count = reports.build_multiverse_cc_qc_report(
                temp_dir,
                atlas_path=atlas_path,
                n_slices=2,
            )

            self.assertEqual(registration_count, 1)
            self.assertEqual(cc_count, 1)
            self.assertTrue(registration_html.is_file())
            self.assertTrue(cc_html.is_file())
            self.assertEqual(
                registration_html.parent,
                Path(temp_dir) / "Report" / "Multiverse" / "Registration",
            )
            self.assertEqual(
                cc_html.parent,
                Path(temp_dir) / "Report" / "Multiverse" / "CC",
            )
            self.assertEqual(
                len(list(registration_html.parent.glob("*.png"))),
                1,
            )
            self.assertEqual(len(list(cc_html.parent.glob("*.png"))), 1)

            registration_text = registration_html.read_text(encoding="utf-8")
            cc_text = cc_html.read_text(encoding="utf-8")
            self.assertIn(mean_path.name, registration_text)
            self.assertIn(atlas_path.name, registration_text)
            self.assertIn("891, 892", cc_text)

    def test_skips_mean_with_mismatched_atlas_affine(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            atlas_affine = np.diag([-0.15, -0.15, 0.15, 1.0])
            atlas_affine[0, 3] = 1.0
            _, atlas_path = self._create_inputs(
                temp_dir,
                atlas_affine=atlas_affine,
            )

            with self.assertLogs(level="WARNING") as captured:
                html_path, count = reports.build_multiverse_registration_qc_report(
                    temp_dir,
                    atlas_path=atlas_path,
                    n_slices=2,
                )

            self.assertIsNone(html_path)
            self.assertEqual(count, 0)
            self.assertIn("different spatial affines", "\n".join(captured.output))


if __name__ == "__main__":
    unittest.main()
