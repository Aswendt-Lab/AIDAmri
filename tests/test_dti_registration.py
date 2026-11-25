"""Exercise atlas resampling and failure handling without NiftyReg binaries."""

import contextlib
import importlib.util
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "bin/2.2_DTIPreProcessing/registration_DTI.py"
spec = importlib.util.spec_from_file_location("registration_DTI", SCRIPT)
registration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(registration)


def option(command, flag):
    return command[command.index(flag) + 1]


class AtlasResamplingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.dti = str(self.folder / "dwi_Bet.nii.gz")
        self.t2 = str(self.folder / "T2_Bet.nii.gz")
        self.template = str(self.folder / "T2_Template.nii.gz")
        self.bspline = str(self.folder / "T2_MatrixBspline.nii")
        self.atlases = [str(self.folder / name) for name in (
            "atlas_split.nii.gz", "atlas_split_parental.nii.gz", "atlas_parental.nii.gz",
        )]
        self.unsplit_atlas = self.folder / "original Allen annotation.nii.gz"
        self.unsplit_atlas.touch()
        (self.folder / "dwi_Bet_mask.nii.gz").touch()
        self.commands = []

    def run_registration(self, fail_on=None):
        def execute(command, **kwargs):
            self.commands.append(command)
            if command[0] == fail_on:
                return subprocess.CompletedProcess(command, 1, "", "synthetic failure")
            if "-res" in command:
                Path(option(command, "-res")).touch()
            return subprocess.CompletedProcess(command, 0, "", "")

        with patch.object(registration.subprocess, "run", side_effect=execute), \
                contextlib.redirect_stdout(io.StringIO()):
            return registration.regABA2DTI(
                self.dti, None, None, self.t2, self.template,
                str(self.folder / "T2_Anno.nii.gz"), *self.atlases,
                self.bspline, str(self.folder),
                allenBrain_anno=str(self.unsplit_atlas),
            )

    def test_atlases_are_sampled_once_from_originals_onto_dti(self):
        output = self.run_registration()
        self.assertEqual(output, str(self.folder / "dwi_Bet_AnnoSplit.nii.gz"))
        composition = self.commands[1]
        self.assertEqual(composition[0], "reg_transform")
        self.assertEqual(option(composition, "-ref"), self.dti)
        self.assertEqual(option(composition, "-ref2"), self.t2)
        first, second, composite = composition[composition.index("-comp") + 1:]
        self.assertEqual(first, option(self.commands[0], "-aff"))
        self.assertEqual(second, self.bspline)

        resamples = [c for c in self.commands if c[0] == "reg_resample"]
        self.assertEqual(len(resamples), 5)  # Four original atlases and the T2 template.
        for atlas in [*self.atlases, str(self.unsplit_atlas)]:
            matches = [c for c in resamples if option(c, "-flo") == atlas]
            self.assertEqual(len(matches), 1)
            command = matches[0]
            self.assertEqual(option(command, "-ref"), self.dti)
            self.assertEqual(option(command, "-trans"), composite)
            self.assertEqual(option(command, "-inter"), "0")
        for command in resamples:
            self.assertNotEqual(option(command, "-flo"), option(command, "-res"))
        template_command = next(c for c in resamples if option(c, "-flo") == self.template)
        self.assertEqual(option(template_command, "-trans"), first)
        anno_command = next(c for c in resamples if option(c, "-flo") == str(self.unsplit_atlas))
        self.assertEqual(option(anno_command, "-res"), str(self.folder / "dwi_Bet_Anno.nii.gz"))

    def test_missing_original_annotation_stops_before_registration(self):
        self.unsplit_atlas.unlink()
        with self.assertRaisesRegex(FileNotFoundError, "Original Allen Brain annotation"):
            self.run_registration()
        self.assertEqual(self.commands, [])

    def test_failed_composition_stops_before_any_resampling(self):
        with self.assertRaisesRegex(RuntimeError, "reg_transform"):
            self.run_registration(fail_on="reg_transform")
        self.assertEqual([c[0] for c in self.commands], ["reg_aladin", "reg_transform"])
        self.assertFalse((self.folder / "DSI_studio").exists())

    def test_failed_rigid_registration_stops_before_composition(self):
        with self.assertRaisesRegex(RuntimeError, "reg_aladin"):
            self.run_registration(fail_on="reg_aladin")
        self.assertEqual([c[0] for c in self.commands], ["reg_aladin"])


if __name__ == "__main__":
    unittest.main()
