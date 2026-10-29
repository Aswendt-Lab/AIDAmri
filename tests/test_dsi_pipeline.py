"""Regression checks for tracking CLI, output selection, and ROI handling."""

import argparse
import contextlib
import io
import os
from pathlib import Path
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin/3.2_DTIConnectivity'))
import dsi_main
import dsi_tools


class TrackingTests(unittest.TestCase):
    def test_named_values_override_mouse_defaults_in_any_order(self):
        first = dsi_main.normalize_track_params(['tract_count=1000', 'step_size=0.05'])
        second = dsi_main.normalize_track_params(['step_size=0.05', 'tract_count=1000'])
        self.assertEqual(first, second)
        self.assertEqual(first['tract_count'], 1000)
        self.assertEqual(first['step_size'], 0.05)
        self.assertEqual(first['max_length'], 15.0)
        self.assertEqual(dsi_main.normalize_track_params(['tract_count=42'])['tract_count'], 42)
        self.assertEqual(dsi_tools.normalize_tracking_params('mouse')['tract_count'], 1000000)

    def test_invalid_cli_values_fail_with_argument_errors(self):
        invalid = [
            ['missing'], ['typo=1'], ['tract_count=1', 'tract_count=2'],
            ['mouse', 'step_size=0.05'], ['1000', '0.05'], ['step_size='],
            ['step_size=NaN'], ['max_length=inf'], ['step_size=0'],
            ['tract_count=1.5'], ['tract_count=-1'], ['check_ending=2'],
            ['smoothing=1.1'], ['turning_angle=181'], ['fa_threshold=-0.1'],
            ['min_length=20', 'max_length=10'],
        ]
        for values in invalid:
            with self.subTest(values=values), self.assertRaises(argparse.ArgumentTypeError):
                dsi_main.normalize_track_params(values)

    def test_all_presets_and_cli_dictionaries_build_tracking_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            fib = Path(tmp) / 'scan.fz'
            fib.touch()
            inputs = [*dsi_tools.TRACKING_PARAM_SETS, 'AIDAmriv2_set',
                      ['tract_count=1000', 'step_size=0.05'],
                      dsi_main.normalize_track_params(['default'])]
            for values in inputs:
                with self.subTest(values=values), patch.object(dsi_tools.subprocess, 'run') as run, \
                        contextlib.redirect_stdout(io.StringIO()):
                    dsi_tools.tracking('dsi_studio', tmp, values, thread_count=3)
                    command = run.call_args.args[0]
                    self.assertIn(f'--source={fib}', command)
                    self.assertIn('--thread_count=3', command)
                    if isinstance(values, list):
                        self.assertIn('--tract_count=1000', command)
                        self.assertIn('--step_size=0.05', command)
                    if values == 'default':
                        self.assertTrue(any(arg.startswith('--parameter_id=') for arg in command))


class PipelineTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.dataset = self.folder / 'dataset'
        self.dataset.mkdir()
        self.image = self.dataset / 'scan.v1_dwi.nii.gz'
        self.affine = np.diag([0.1, 0.1, 0.5, 1.0])
        nib.save(nib.Nifti1Image(np.zeros((4, 4, 2, 3)), self.affine), self.image)
        self.mask = self.dataset / 'mask.nii.gz'
        self.mask.touch()
        self.btable = self.dataset / 'btable.txt'
        self.btable.touch()
        self.cwd = os.getcwd()
        self.addCleanup(os.chdir, self.cwd)
        stdout = contextlib.redirect_stdout(io.StringIO())
        stdout.__enter__()
        self.addCleanup(stdout.__exit__, None, None, None)

    def test_relative_srcgen_uses_current_fib_and_keeps_cwd(self):
        fib_dir = self.dataset / 'fib_map'
        fib_dir.mkdir()
        old = fib_dir / 'aaa_stale.fz'
        old.write_text('old reconstruction')
        commands = []

        def run(command, **kwargs):
            commands.append(command)
            self.assertEqual(kwargs['cwd'], str(self.dataset))
            output = next((arg.split('=', 1)[1] for arg in command if arg.startswith('--output=')), None)
            if output:
                self.assertTrue(Path(output).is_absolute())
                Path(output).touch()
            if '--action=exp' in command:
                source = next(arg.split('=', 1)[1] for arg in command if arg.startswith('--source='))
                metric = next(arg.split('=', 1)[1] for arg in command if arg.startswith('--export='))
                nib.save(nib.Nifti1Image(np.zeros((4, 4, 2)), self.affine), source + f'.{metric}.nii.gz')

        os.chdir(self.folder)
        with patch.object(dsi_tools.time, 'sleep'), patch.object(dsi_tools.subprocess, 'run', side_effect=run):
            voxel = dsi_tools.srcgen('./dsi_studio', 'dataset/' + self.image.name,
                                    'dataset/mask.nii.gz', 'dataset/' + self.image.name, 'btable.txt')
        self.assertEqual(os.getcwd(), str(self.folder))
        self.assertAlmostEqual(voxel, 0.1)
        exports = [command for command in commands if '--action=exp' in command]
        for command in exports:
            self.assertIn(f'--source={fib_dir}/scan.v1_dwi.fz', command)
            self.assertEqual(command[0], str(self.folder / 'dsi_studio'))
        self.assertEqual(old.read_text(), 'old reconstruction')
        metrics = list((self.dataset / 'DSI_studio').glob('*.nii.gz'))
        self.assertEqual(len(metrics), 4)
        self.assertEqual(nib.aff2axcodes(nib.load(metrics[0]).affine), ('L', 'I', 'P'))

    def test_srcgen_failure_preserves_cwd(self):
        os.chdir(self.folder)
        with patch.object(dsi_tools.time, 'sleep'), \
                patch.object(dsi_tools.subprocess, 'run', side_effect=RuntimeError('DSI failed')), \
                self.assertRaisesRegex(RuntimeError, 'DSI failed'):
            dsi_tools.srcgen('dsi_studio', 'dataset/' + self.image.name, 'mask.nii.gz',
                             'dataset/' + self.image.name, 'btable.txt')
        self.assertEqual(os.getcwd(), str(self.folder))

    def test_legacy_fallback_rejects_stale_and_ambiguous_files(self):
        directory = str(self.dataset)
        stale = self.dataset / 'scan.src.gz.dti.fib.gz'
        stale.touch()
        before = dsi_tools.snapshot_dsi_outputs([directory], 'scan', ('.fib.gz',))
        expected = str(self.dataset / 'scan.fib.gz')
        with self.assertRaises(FileNotFoundError):
            dsi_tools.select_dsi_output(expected, before, [directory], 'scan', ('.fib.gz',))
        current = self.dataset / 'scan.src.gz.gqi.fib.gz'
        current.touch()
        self.assertEqual(dsi_tools.select_dsi_output(expected, before, [directory], 'scan', ('.fib.gz',)), str(current))
        (self.dataset / 'scan.other.fib.gz').touch()
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            dsi_tools.select_dsi_output(expected, before, [directory], 'scan', ('.fib.gz',))

    def test_unchanged_expected_output_does_not_hide_new_legacy_output(self):
        directory = str(self.dataset)
        expected = self.dataset / 'scan.fz'
        expected.write_text('stale')
        before = dsi_tools.snapshot_dsi_outputs([directory], 'scan', ('.fz',))
        current = self.dataset / 'scan.sz.dti.fz'
        current.write_text('current')
        self.assertEqual(dsi_tools.select_dsi_output(str(expected), before, [directory], 'scan', ('.fz',)), str(current))

    def test_relative_gradients_are_anchored_for_dsi(self):
        bval = self.dataset / 'scan.bval'
        bvec = self.dataset / 'scan.bvec'
        bval.touch()
        bvec.touch()
        os.chdir(self.folder)
        with patch.object(dsi_tools.time, 'sleep'), \
                patch.object(dsi_tools.subprocess, 'run', side_effect=RuntimeError('stop after command')) as run, \
                self.assertRaises(RuntimeError):
            dsi_tools.srcgen('dsi_studio', 'dataset/' + self.image.name, 'mask.nii.gz',
                             'dataset/' + self.image.name, None,
                             gradient_pair={'bval': 'dataset/scan.bval', 'bvec': 'scan.bvec'})
        self.assertIn(f'--bval={bval}', run.call_args.args[0])
        self.assertIn(f'--bvec={bvec}', run.call_args.args[0])

    def test_resampling_preserves_stroke_and_atlas_label_files(self):
        fib = self.dataset / 'fib_map'
        fib.mkdir()
        (fib / 'scan.fz').touch()
        (fib / 'scan.fz.trk.gz').touch()
        for stem in ('scan_Stroke_mask_anno', 'scan_AnnoSplit_parental'):
            roi = self.dataset / (stem + '.nii.gz')
            nib.save(nib.Nifti1Image(np.ones((4, 4, 2)), self.affine), roi)
        (self.dataset / 'scan_AnnoSplit_parental.txt').write_text('1 Motor\n')
        os.chdir(self.folder)
        for stem, expected_label in [('scan_Stroke_mask_anno', '1 Stroke_Motor\n'),
                                     ('scan_AnnoSplit_parental', '1 Motor\n')]:
            with self.subTest(stem=stem), patch.object(dsi_tools.time, 'sleep'), \
                    patch.object(dsi_tools.subprocess, 'run') as run:
                dsi_tools.connectivity('dsi_studio', 'dataset/fib_map', f'dataset/{stem}.nii.gz',
                                       'dataset', 'connectivity', make_isotropic=0.1)
                roi = self.dataset / 'connectivity' / (stem + '_resampled.nii.gz')
                label = self.dataset / 'connectivity' / (stem + '_resampled.txt')
                self.assertEqual(label.read_text(), expected_label)
                analysis = [call.args[0] for call in run.call_args_list if '--action=ana' in call.args[0]]
                self.assertEqual(len(analysis), 2)
                self.assertIn(f'--connectivity={roi}', analysis[0])

    def test_motion_correction_failure_cleans_only_its_temporary_files(self):
        par = self.dataset / 'mcf_Folder'
        par.mkdir()
        old_temp = self.dataset / 'temp'
        old_temp.mkdir()
        marker = old_temp / 'keep.txt'
        marker.write_text('keep')
        def fail(**kwargs):
            raise RuntimeError('FSL failed')

        split = SimpleNamespace(run=fail)
        with patch.object(dsi_tools.fsl, 'Split', return_value=split), self.assertRaisesRegex(RuntimeError, 'FSL failed'):
            dsi_tools.fsl_SeparateSliceMoCo(str(self.image), str(par))
        self.assertEqual(os.getcwd(), self.cwd)
        self.assertEqual(marker.read_text(), 'keep')
        self.assertFalse((self.dataset / 'fslScaleTemp.nii.gz').exists())
        self.assertEqual(list(self.dataset.glob('dsi-moco-*')), [])

    def test_motion_correction_merges_only_current_slices_and_restores_scale(self):
        par = self.dataset / 'mcf_Folder'
        par.mkdir()
        unrelated = par / 'scan.v1_dwi_old.nii.gz'
        unrelated.write_text('keep')
        scaled = self.dataset / 'fslScaleTemp.nii.gz'

        def split_run(cwd):
            for index in range(2):
                (Path(cwd) / f'scan.v1_dwi{index:04d}.nii.gz').touch()

        def motion(**kwargs):
            def run(cwd):
                shutil.copyfile(kwargs['in_file'], kwargs['out_file'])
            return SimpleNamespace(run=run)

        def merge(**kwargs):
            self.assertEqual(len(kwargs['in_files']), 2)
            self.assertNotIn(str(unrelated), kwargs['in_files'])
            return SimpleNamespace(run=lambda cwd: shutil.copyfile(scaled, kwargs['merged_file']))

        with patch.object(dsi_tools.fsl, 'Split', return_value=SimpleNamespace(run=split_run)), \
                patch.object(dsi_tools.fsl.preprocess, 'MCFLIRT', side_effect=motion), \
                patch.object(dsi_tools.fsl, 'Merge', side_effect=merge):
            output = dsi_tools.fsl_SeparateSliceMoCo(str(self.image), str(par))
        self.assertEqual(Path(output).name, 'scan.v1_dwi_mcf.nii.gz')
        np.testing.assert_allclose(nib.load(output).affine, self.affine)
        self.assertEqual(unrelated.read_text(), 'keep')
        self.assertEqual(os.getcwd(), self.cwd)


if __name__ == '__main__':
    unittest.main()
