"""Check denoiser integration and gradient lookup with small synthetic images."""

import contextlib
import io
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin/2.2_DTIPreProcessing"))
sys.path.insert(0, str(ROOT / "bin"))
import preProcessing_DTI as preprocessing
import averageb0
import batchProc


class DenoisingTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory(prefix="dwi denoising ")
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)
        self.input = self.folder / "scan_dwi.nii.gz"
        self.data = np.zeros((5, 5, 5, 4), dtype=np.float32)
        self.data[..., 0] = 2
        self.data[..., 1] = 6
        self.data[..., 2:] = 20
        self.affine = np.diag([-0.2, -0.3, 0.4, 1])
        nib.save(nib.Nifti1Image(self.data, self.affine), self.input)
        self.bvals = self.folder / "scan_dwi.bval"
        np.savetxt(self.bvals, [0, 0, 1000, 1000])
        self.addCleanup(os.chdir, os.getcwd())

    def test_mrtrix_output_geometry_and_noise_name(self):
        def run(command, **kwargs):
            self.assertTrue(kwargs["check"])
            self.assertEqual(command[1], str(self.input))
            self.assertEqual(Path(command[2]).name, "scan_dwiMR3Denoised.nii.gz")
            noise = command[command.index("-noise") + 1]
            self.assertEqual(Path(noise).name, "scan_dwiMR3Noise.nii.gz")
            nib.save(nib.load(command[1]), command[2])
            nib.save(nib.Nifti1Image(self.data[..., 0], self.affine), noise)
            return subprocess.CompletedProcess(command, 0, "done\n")

        with patch.object(preprocessing.shutil, "which", return_value="/mock/dwidenoise"), \
                patch.object(preprocessing.subprocess, "run", side_effect=run), \
                contextlib.redirect_stdout(io.StringIO()):
            output = preprocessing.denoise_mrtrix3(str(self.input), str(self.folder))
        self.assertEqual(nib.load(output).shape, self.data.shape)
        np.testing.assert_allclose(nib.load(output).affine, self.affine)
        np.testing.assert_array_equal(nib.load(self.input).get_fdata(), self.data)

    def test_missing_mrtrix_and_command_failure_are_not_silently_ignored(self):
        with patch.object(preprocessing.shutil, "which", return_value=None), \
                self.assertRaisesRegex(FileNotFoundError, "MRtrix3"):
            preprocessing.denoise_mrtrix3(str(self.input), str(self.folder))
        failure = subprocess.CalledProcessError(1, ["dwidenoise"], output="failed")
        with patch.object(preprocessing.shutil, "which", return_value="dwidenoise"), \
                patch.object(preprocessing.subprocess, "run", side_effect=failure), \
                contextlib.redirect_stdout(io.StringIO()), \
                self.assertRaises(subprocess.CalledProcessError):
            preprocessing.denoise_mrtrix3(str(self.input), str(self.folder))

    def test_mrtrix_rejects_3d_input_before_running(self):
        nib.save(nib.Nifti1Image(self.data[..., 0], self.affine), self.input)
        with patch.object(preprocessing.shutil, "which", return_value="dwidenoise"), \
                patch.object(preprocessing.subprocess, "run") as run, \
                self.assertRaisesRegex(ValueError, "4D"):
            preprocessing.denoise_mrtrix3(str(self.input), str(self.folder))
        run.assert_not_called()

    def test_patch2self_uses_short_name_and_retains_volume_order(self):
        dipy = types.ModuleType("dipy")
        dipy.denoise = types.ModuleType("dipy.denoise")
        module = types.ModuleType("dipy.denoise.patch2self")
        dipy.denoise.patch2self = module
        module.patch2self = Mock(return_value=self.data)
        with patch.dict(sys.modules, {
            "dipy": dipy, "dipy.denoise": dipy.denoise, "dipy.denoise.patch2self": module,
        }), patch.object(preprocessing.subprocess, 'run') as noise_run:
            output = preprocessing.denoise_patch2self(str(self.input), str(self.folder))
        self.assertEqual(Path(output).name, "scan_dwiP2SDenoised.nii.gz")
        np.testing.assert_array_equal(module.patch2self.call_args.args[1], [0, 0, 1000, 1000])
        np.testing.assert_array_equal(nib.load(output).get_fdata(), self.data)

        noise_run.assert_not_called()
        noise = nib.load(self.folder / 'scan_dwiP2SNoise.nii.gz')
        np.testing.assert_array_equal(noise.get_fdata(), 0)

    def test_patch2self_noise_is_signed_4d_difference_without_mrtrix(self):
        expected = np.zeros_like(self.data)
        expected[0, ...] = [0, 2, -3, 4]
        expected[1, ...] = [-1, 0, 5, -6]
        denoised = self.data - expected
        source = nib.load(self.input)
        source.header.set_xyzt_units('mm', 'sec')
        nib.save(source, self.input)
        output_folder = self.folder / 'outputs'
        output_folder.mkdir()
        dipy = types.ModuleType('dipy')
        dipy.denoise = types.ModuleType('dipy.denoise')
        module = types.ModuleType('dipy.denoise.patch2self')
        dipy.denoise.patch2self = module
        module.patch2self = Mock(return_value=denoised)
        with patch.dict(sys.modules, {
            'dipy': dipy, 'dipy.denoise': dipy.denoise, 'dipy.denoise.patch2self': module,
        }), patch.object(preprocessing.shutil, 'which', return_value=None), \
                patch.object(preprocessing.subprocess, 'run') as run:
            output = preprocessing.denoise_patch2self(str(self.input), str(output_folder))
        run.assert_not_called()
        noise = nib.load(output_folder / 'scan_dwiP2SNoise.nii.gz')
        self.assertEqual(noise.shape, self.data.shape)
        self.assertEqual(noise.get_data_dtype(), np.dtype(np.float32))
        self.assertEqual(noise.header.get_xyzt_units(), ('mm', 'sec'))
        self.assertIn(b'original minus denoised', noise.header['descrip'].tobytes())
        np.testing.assert_array_equal(noise.get_fdata(), expected)
        np.testing.assert_allclose(noise.affine, self.affine)
        for transform in (noise.get_qform, noise.get_sform):
            affine, code = transform(coded=True)
            self.assertEqual(code, 1)
            np.testing.assert_allclose(affine, self.affine)
        np.testing.assert_array_equal(nib.load(self.input).get_fdata(), self.data)
        np.testing.assert_array_equal(nib.load(output).get_fdata(), denoised)
        np.testing.assert_array_equal(nib.load(output).get_fdata() + noise.get_fdata(), self.data)

    def test_b0_averaging_uses_original_gradients_for_all_suffixes(self):
        for extension in (".bval", ".btable"):
            if extension == ".btable":
                self.bvals.unlink()
                np.savetxt(self.folder / "scan_dwi.btable", [[0, 0, 1000, 1000]] * 4)
            for suffix in ("P2SDenoised", "MR3Denoised", "Patch2SelfDenoised"):
                with self.subTest(extension=extension, suffix=suffix):
                    if extension == ".btable":
                        self.bvals.unlink(missing_ok=True)
                    denoised = self.folder / f"scan_dwi{suffix}.nii.gz"
                    nib.save(nib.load(self.input), denoised)
                    output = averageb0.averageb0(str(denoised))
                    np.testing.assert_array_equal(nib.load(output).get_fdata(), 4)
                    np.testing.assert_allclose(nib.load(output).affine, self.affine)

    def test_batch_passes_selected_output_and_supports_legacy_patch2self(self):
        for method, suffix in (("patch2self", "P2SDenoised"),
                               ("mrtrix3", "MR3Denoised"),
                               ("patch2self", "Patch2SelfDenoised")):
            with self.subTest(method=method, suffix=suffix):
                output = self.folder / f"scan_dwi{suffix}.nii.gz"
                output.touch()
                (self.folder / 'scan_dwiP2SNoise.nii.gz').touch()
                with patch.object(batchProc, "run_subprocess", return_value=0) as run:
                    batchProc.executeScripts(self.folder, "dwi", "process", {"dwi_denoiser": method})
                command = shlex.split(run.call_args.args[0])
                self.assertEqual(command[command.index("-i") + 1], str(output))
                output.unlink()

    def test_batch_forwards_mrtrix_option_to_preprocessing(self):
        with patch.object(batchProc, "run_subprocess", return_value=0) as run:
            batchProc.executeScripts(self.folder, "dwi", "preprocess", {"dwi_denoiser": "mrtrix3"})
        self.assertIn("--denoiser mrtrix3", run.call_args.args[0])


if __name__ == "__main__":
    unittest.main()
