"""Exercise batch command construction through dsi_main and real DSI helpers."""

import atexit
import contextlib
import io
import os
from pathlib import Path
import runpy
import shlex
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np
from scipy.io import loadmat, savemat

SCRIPT_DIR = Path(__file__).resolve().parents[1] / 'bin'
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPT_DIR / '3.2_DTIConnectivity'))
import batchProc
import dsi_tools
from common import artifact_manifest


class DsiBatchIntegrationTests(unittest.TestCase):
    def test_batch_calls_main_and_tracking_for_presets_custom_values_and_legacy(self):
        scenarios = [
            (None, '0', False, False, 'dti'),
            (['mouse'], '0', False, False, 'dti'),
            (['AIDAmriv2_set'], '0.2', False, False, 'gqi'),
            (['tract_count=1000', 'step_size=0.05'], 'auto', True, False, 'dti'),
            (('step_size=0.05', 'tract_count=1000'), 'auto', False, True, 'gqi'),
        ]
        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        for values, isotropic, legacy, denoised, method in scenarios:
            with self.subTest(values=values, isotropic=isotropic, legacy=legacy, denoised=denoised), \
                    tempfile.TemporaryDirectory(prefix='dsi batch ') as tmp:
                folder = Path(tmp) / 'sub-01/ses-01/dwi'
                folder.mkdir(parents=True)
                affine = np.diag([0.1, 0.1, 0.5, 1.0])
                image = folder / 'scan_dwi.nii.gz'
                nib.save(nib.Nifti1Image(np.zeros((4, 4, 2, 3)), affine), image)
                (folder / 'scan_dwi.bval').touch()
                (folder / 'scan_dwi.bvec').touch()
                (folder / 'scanBet_mask.nii.gz').touch()
                seed = folder / 'scan_AnnoSplit_parental.nii.gz'
                nib.save(nib.Nifti1Image(np.ones((4, 4, 2)), affine), seed)
                (folder / 'scan_AnnoSplit_parental.txt').write_text('1 Motor\n')
                if denoised:
                    shutil.copyfile(image, folder / 'scan_P2SDenoised.nii.gz')
                calls = []
                source_inputs = []

                def reconstruct(dsi_studio, dir_in, dir_msk, dir_out, b_table,
                                recon_method='dti', vivo='in_vivo', make_isotropic=0,
                                legacy=False, gradient_pair=None):
                    source_inputs.append(dir_in)
                    self.assertEqual(recon_method, method)
                    self.assertEqual(gradient_pair['bval'], str(folder / 'scan_dwi.bval'))
                    fib_dir = folder / 'fib_map'
                    fib_dir.mkdir()
                    (fib_dir / ('scan.fib.gz' if legacy else 'scan.fz')).touch()
                    (folder / 'DSI_studio').mkdir()
                    return 0.1

                def external(command, **kwargs):
                    calls.append(command)
                    if '--action=trk' in command:
                        output = next(arg.split('=', 1)[1] for arg in command if arg.startswith('--output='))
                        Path(output).touch()
                    elif '--action=ana' in command:
                        tract = next(arg.split('=', 1)[1] for arg in command if arg.startswith('--tract='))
                        roi = next(arg.split('=', 1)[1] for arg in command if arg.startswith('--connectivity='))
                        savemat(tract + '.' + dsi_tools.strip_nifti_suffix(roi) + '.connectivity.mat', {
                            'name': np.frombuffer(b'Motor\n\0', dtype=np.uint8)[None, :],
                            'number of tracts t2r': np.array([[1000]], dtype=np.float32),
                        }, format='4')
                    elif command[0] == 'flirt':
                        shutil.copyfile(command[command.index('-in') + 1], command[command.index('-out') + 1])

                def child(command, datatype, step):
                    argv = shlex.split(command)
                    self.assertEqual(argv[:2], ['python', 'dsi_main.py'])
                    with patch.object(sys, 'argv', argv[1:]):
                        runpy.run_path(str(SCRIPT_DIR / '3.2_DTIConnectivity/dsi_main.py'), run_name='__main__')
                    return 0

                cfg = {
                    'num_processes': 3, 'dsi_track_param': values,
                    'dsi_make_isotropic': isotropic, 'dsi_legacy': legacy,
                    'dsi_skip_motion_correction': True, 'dsi_recon_method': method,
                    'dwi_denoiser': 'patch2self' if denoised else None,
                }
                with patch.object(batchProc, 'run_subprocess', side_effect=child), \
                        patch.object(dsi_tools, 'srcgen', autospec=True, side_effect=reconstruct), \
                        patch.object(dsi_tools.subprocess, 'run', side_effect=external), \
                        patch.object(artifact_manifest, 'start_output_tracking'), \
                        patch.object(atexit, 'register'), \
                        patch.dict(os.environ, {'AIDAMRI_DISABLE_SCRIPT_LOG': '1'}), \
                        contextlib.redirect_stdout(io.StringIO()):
                    result = batchProc.executeScripts(folder, 'dwi', 'process', cfg)
                self.assertEqual(result, 0)
                tracking = next(command for command in calls if '--action=trk' in command)
                self.assertIn('--thread_count=3', tracking)
                self.assertTrue(any(arg.endswith('scan.fib.gz' if legacy else 'scan.fz') for arg in tracking))
                if isinstance(values, (tuple, list)) and '=' in values[0]:
                    self.assertIn('--tract_count=1000', tracking)
                    self.assertIn('--step_size=0.05', tracking)
                if denoised:
                    self.assertTrue(source_inputs[0].endswith('scan_P2SDenoised.nii.gz'))
                if isotropic == 'auto':
                    resample = next(command for command in calls if command[0] == 'flirt')
                    self.assertEqual(float(resample[resample.index('-applyisoxfm') + 1]), 0.1)
                results = list((folder / 'connectivity').glob('*.mat'))
                self.assertEqual(len(results), 2)
                self.assertTrue(any('.all.pass.connectivity.mat' in path.name for path in results))
                self.assertTrue(any('.all.end.connectivity.mat' in path.name for path in results))
                self.assertIn('region_names', loadmat(results[0]))
                self.assertEqual(list(folder.rglob('*.csv')), [])
                os.chdir(cwd)


if __name__ == '__main__':
    unittest.main()
