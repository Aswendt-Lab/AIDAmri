"""Tests for transforming processed 4D fMRI data into SIGMA space."""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import nibabel as nib
import numpy as np


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
MODULE_PATH = BIN_DIR / "Create_multiverse_output.py"
sys.path.insert(0, str(BIN_DIR))


def load_module():
    spec = importlib.util.spec_from_file_location(
        "create_multiverse_output_transform_test",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save_nifti(path, shape, affine=None):
    if affine is None:
        affine = np.eye(4)
    nib.save(nib.Nifti1Image(np.zeros(shape, dtype=np.float32), affine), path)


class MultiverseTransformTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_masks_refines_and_applies_corrected_composite_transform(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            func_root = Path(temp_dir) / "func"
            func_folder = func_root / "rs-fMRI_niiData"
            func_folder.mkdir(parents=True)

            prefix = "sub-01_task-rest_bold_EPIMPSmoothBet"
            reference = func_root / f"{prefix}.nii.gz"
            bet_mask = func_root / f"{prefix}_mask.nii.gz"
            composite = func_root / f"{prefix}Matrixcomp_rsfMRI.nii.gz"
            template = Path(temp_dir) / "SIGMA_template.nii.gz"
            atlas = Path(temp_dir) / self.module.SIGMA_ATLAS_FILENAME
            input_file = func_folder / "sub-01_task-rest_bold_EPI_mcf_f.nii.gz"

            save_nifti(reference, (4, 5, 6))
            mask_data = np.zeros((4, 5, 6), dtype=np.uint8)
            mask_data[1:3, 1:4, 1:5] = 1
            nib.save(nib.Nifti1Image(mask_data, np.eye(4)), bet_mask)
            sigma_affine = np.diag([0.15, 0.15, 0.15, 1.0])
            save_nifti(template, (8, 10, 12), sigma_affine)
            save_nifti(atlas, (8, 10, 12), sigma_affine)
            input_data = np.ones((4, 5, 6, 3), dtype=np.float32)
            nib.save(nib.Nifti1Image(input_data, np.eye(4)), input_file)
            composite.write_bytes(b"composite")

            def fake_niftyreg_run(command, check):
                self.assertTrue(check)
                output_path = Path(command[-1])
                if command[0] == "reg_transform":
                    output_path.write_bytes(b"transform")
                elif command[0] == "reg_resample":
                    floating = nib.load(command[command.index("-flo") + 1])
                    reference_img = nib.load(command[command.index("-ref") + 1])
                    output_path = Path(command[command.index("-res") + 1])
                    output_shape = reference_img.shape[:3] + floating.shape[3:]
                    save_nifti(output_path, output_shape, reference_img.affine)
                elif command[0] == "reg_aladin":
                    result_path = Path(command[command.index("-res") + 1])
                    affine_path = Path(command[command.index("-aff") + 1])
                    reference_img = nib.load(command[command.index("-ref") + 1])
                    save_nifti(result_path, reference_img.shape, reference_img.affine)
                    affine_path.write_text("identity", encoding="utf-8")
                return mock.DEFAULT

            with mock.patch.object(
                self.module.subprocess,
                "run",
                side_effect=fake_niftyreg_run,
            ) as run:
                outputs = self.module.apply_inverse_composite_transformation(
                    files_list=[str(input_file)],
                    func_folder=str(func_folder),
                    sigma_template_address=str(template),
                )

            output_dir = func_folder / self.module.MULTIVERSE_OUTPUT_FOLDER
            inverse = output_dir / f"{prefix}Matrixcomp_rsfMRI_inv.nii.gz"
            provisional = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template.nii.gz"
            )
            provisional_mean = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_temporal_mean.nii.gz"
            )
            second_registration = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_2nd_reg.nii.gz"
            )
            second_affine = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_2nd_aff.txt"
            )
            corrected_transform = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_corrected_trans.nii.gz"
            )
            corrected = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_corrected.nii.gz"
            )
            corrected_mean = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_temporal_mean_corrected.nii.gz"
            )
            corrected_resampled = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_corrected_resampled_0p3mm.nii.gz"
            )
            corrected_mean_resampled = output_dir / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_temporal_mean_corrected_"
                "resampled_0p3mm.nii.gz"
            )
            masked_input = output_dir / "sub-01_task-rest_bold_EPI_mcf_f_BET.nii.gz"
            self.assertEqual(outputs, [str(corrected_resampled)])
            masked_data = nib.load(masked_input).get_fdata()
            self.assertEqual(masked_data.shape, (4, 5, 6, 3))
            self.assertTrue(np.all(masked_data[mask_data == 0] == 0))
            self.assertTrue(np.all(masked_data[mask_data > 0] == 1))
            self.assertEqual(nib.load(provisional).shape, (8, 10, 12, 3))
            self.assertEqual(nib.load(provisional_mean).shape, (8, 10, 12))
            self.assertEqual(nib.load(second_registration).shape, (8, 10, 12))
            self.assertTrue(second_affine.is_file())
            self.assertTrue(corrected_transform.is_file())
            self.assertEqual(nib.load(corrected).shape, (8, 10, 12, 3))
            self.assertEqual(nib.load(corrected_mean).shape, (8, 10, 12))
            self.assertEqual(
                nib.load(corrected_resampled).header.get_zooms()[:3],
                (0.3, 0.3, 0.3),
            )
            self.assertEqual(
                nib.load(corrected_mean_resampled).header.get_zooms()[:3],
                (0.3, 0.3, 0.3),
            )
            self.assertEqual(
                nib.load(corrected_resampled).shape,
                (4, 5, 6, 3),
            )
            self.assertEqual(
                nib.load(corrected_mean_resampled).shape,
                (4, 5, 6),
            )
            expected_resampled_affine = sigma_affine.copy()
            expected_resampled_affine[:3, :3] *= 2
            expected_resampled_affine[:3, 3] = 0.075
            np.testing.assert_allclose(
                nib.load(corrected_resampled).affine,
                expected_resampled_affine,
            )
            self.assertEqual(
                run.call_args_list,
                [
                    mock.call(
                        [
                            "reg_transform",
                            "-invNrr",
                            str(composite),
                            str(template),
                            str(inverse),
                        ],
                        check=True,
                    ),
                    mock.call(
                        [
                            "reg_resample",
                            "-ref",
                            str(template),
                            "-flo",
                            str(masked_input),
                            "-trans",
                            str(inverse),
                            "-res",
                            str(provisional),
                            "-inter",
                            "3",
                        ],
                        check=True,
                    ),
                    mock.call(
                        [
                            "reg_aladin",
                            "-ref",
                            str(template),
                            "-flo",
                            str(provisional_mean),
                            "-rigOnly",
                            "-res",
                            str(second_registration),
                            "-aff",
                            str(second_affine),
                        ],
                        check=True,
                    ),
                    mock.call(
                        [
                            "reg_transform",
                            "-ref",
                            str(template),
                            "-comp",
                            str(second_affine),
                            str(inverse),
                            str(corrected_transform),
                        ],
                        check=True,
                    ),
                    mock.call(
                        [
                            "reg_resample",
                            "-ref",
                            str(template),
                            "-flo",
                            str(masked_input),
                            "-trans",
                            str(corrected_transform),
                            "-res",
                            str(corrected),
                            "-inter",
                            "3",
                        ],
                        check=True,
                    ),
                ],
            )

    def test_rejects_non_0p15mm_grid_before_resampling(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir) / "corrected.nii.gz"
            output_file = Path(temp_dir) / "resampled.nii.gz"
            wrong_affine = np.diag([0.15, 0.2, 0.15, 1.0])
            save_nifti(input_file, (8, 10, 12), wrong_affine)

            with self.assertRaisesRegex(
                ValueError,
                "Expected a 0.15 mm isotropic input grid",
            ):
                self.module.downsample_nifti_by_block_mean(
                    str(input_file),
                    str(output_file),
                )

            self.assertFalse(output_file.exists())

    def test_downsamples_each_4d_volume_by_2x2x2_block_mean(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir) / "corrected.nii.gz"
            output_file = Path(temp_dir) / "resampled.nii.gz"
            affine = np.diag([0.15, 0.15, 0.15, 1.0])
            x, y, z = np.indices((4, 4, 4))
            first_volume = 100 * x + 10 * y + z
            input_data = np.stack(
                [first_volume, first_volume + 1000],
                axis=3,
            ).astype(np.float32)
            nib.save(nib.Nifti1Image(input_data, affine), input_file)

            self.module.downsample_nifti_by_block_mean(
                str(input_file),
                str(output_file),
            )

            output_img = nib.load(output_file)
            output_data = output_img.get_fdata(dtype=np.float32)
            self.assertEqual(output_data.shape, (2, 2, 2, 2))
            self.assertEqual(output_data[0, 0, 0, 0], 55.5)
            self.assertEqual(output_data[1, 1, 1, 0], 277.5)
            self.assertEqual(output_data[0, 0, 0, 1], 1055.5)
            self.assertEqual(output_data[1, 1, 1, 1], 1277.5)

    def test_rejects_fmri_on_a_different_spatial_grid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            func_root = Path(temp_dir) / "func"
            func_folder = func_root / "rs-fMRI_niiData"
            func_folder.mkdir(parents=True)

            prefix = "sub-01_task-rest_bold_EPIMPSmoothBet"
            reference = func_root / f"{prefix}.nii.gz"
            bet_mask = func_root / f"{prefix}_mask.nii.gz"
            composite = func_root / f"{prefix}Matrixcomp_rsfMRI.nii.gz"
            template = Path(temp_dir) / "SIGMA_template.nii.gz"
            input_file = func_folder / "sub-01_task-rest_bold_EPI_mcf_f.nii.gz"

            save_nifti(reference, (4, 5, 6))
            save_nifti(bet_mask, (4, 5, 6))
            save_nifti(template, (7, 8, 9))
            save_nifti(input_file, (4, 5, 7, 3))
            composite.write_bytes(b"composite")

            with (
                mock.patch.object(self.module.subprocess, "run") as run,
                self.assertRaisesRegex(ValueError, "different spatial shapes"),
            ):
                self.module.apply_inverse_composite_transformation(
                    files_list=[str(input_file)],
                    func_folder=str(func_folder),
                    sigma_template_address=str(template),
                )

            run.assert_not_called()

    def test_requires_companion_bet_mask(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            func_root = Path(temp_dir) / "func"
            func_folder = func_root / "rs-fMRI_niiData"
            func_folder.mkdir(parents=True)

            prefix = "sub-01_task-rest_bold_EPIMPSmoothBet"
            reference = func_root / f"{prefix}.nii.gz"
            composite = func_root / f"{prefix}Matrixcomp_rsfMRI.nii.gz"
            template = Path(temp_dir) / "SIGMA_template.nii.gz"
            input_file = func_folder / "sub-01_task-rest_bold_EPI_mcf_f.nii.gz"

            save_nifti(reference, (4, 5, 6))
            save_nifti(template, (7, 8, 9))
            save_nifti(input_file, (4, 5, 6, 3))
            composite.write_bytes(b"composite")

            with (
                mock.patch.object(self.module.subprocess, "run") as run,
                self.assertRaisesRegex(FileNotFoundError, "BET mask"),
            ):
                self.module.apply_inverse_composite_transformation(
                    files_list=[str(input_file)],
                    func_folder=str(func_folder),
                    sigma_template_address=str(template),
                )

            run.assert_not_called()

    def test_accepts_minor_ants_affine_normalisation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir) / "fmri.nii.gz"
            reference = Path(temp_dir) / "reference.nii.gz"
            affine = np.diag([-0.4, -0.4, -1.1, 1.0])
            reference_affine = affine.copy()
            reference_affine[1, 2] = 0.0006

            save_nifti(input_file, (64, 64, 18, 3), affine)
            save_nifti(reference, (64, 64, 18), reference_affine)

            with mock.patch("builtins.print") as print_mock:
                self.module.validate_spatial_geometry(
                    str(input_file),
                    str(reference),
                )

            self.assertIn(
                "Accepting a minor fMRI/reference affine difference",
                print_mock.call_args.args[0],
            )

    def test_rejects_material_affine_difference(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir) / "fmri.nii.gz"
            reference = Path(temp_dir) / "reference.nii.gz"
            affine = np.diag([-0.4, -0.4, -1.1, 1.0])
            shifted_affine = affine.copy()
            shifted_affine[0, 3] = 0.2

            save_nifti(input_file, (64, 64, 18, 3), affine)
            save_nifti(reference, (64, 64, 18), shifted_affine)

            with self.assertRaisesRegex(
                ValueError,
                "materially different spatial affines",
            ):
                self.module.validate_spatial_geometry(
                    str(input_file),
                    str(reference),
                )

    def test_rejects_sigma_template_with_wrong_atlas_affine(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            template = Path(temp_dir) / "SIGMA_InVivo_Brain_Template_Masked.nii"
            atlas = Path(temp_dir) / self.module.SIGMA_ATLAS_FILENAME
            template_affine = np.eye(4)
            atlas_affine = np.eye(4)
            atlas_affine[1, 3] = 11.045

            save_nifti(template, (128, 128, 218), template_affine)
            save_nifti(atlas, (128, 128, 218), atlas_affine)

            with self.assertRaisesRegex(
                ValueError,
                "does not use the geometry of the SIGMA atlas",
            ):
                self.module.validate_sigma_template_geometry(str(template))

    def test_accepts_sigma_template_matching_atlas_geometry(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            template = Path(temp_dir) / "SIGMA_InVivo_Brain_Template_Masked.nii.gz"
            atlas = Path(temp_dir) / self.module.SIGMA_ATLAS_FILENAME
            affine = np.diag([-0.15, -0.15, -0.15, 1.0])

            save_nifti(template, (128, 128, 218), affine)
            save_nifti(atlas, (128, 128, 218), affine)

            self.module.validate_sigma_template_geometry(str(template))

    def test_rejects_ambiguous_composite_transforms(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            func_root = Path(temp_dir) / "func"
            func_folder = func_root / "rs-fMRI_niiData"
            func_folder.mkdir(parents=True)
            (func_root / "firstMatrixcomp_rsfMRI.nii.gz").write_bytes(b"first")
            (func_root / "secondMatrixcomp_rsfMRI.nii.gz").write_bytes(b"second")

            with self.assertRaisesRegex(RuntimeError, "found 2"):
                self.module.apply_inverse_composite_transformation(
                    files_list=[],
                    func_folder=str(func_folder),
                    sigma_template_address="template.nii.gz",
                )


if __name__ == "__main__":
    unittest.main()
