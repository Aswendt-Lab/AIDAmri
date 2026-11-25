"""Check final conversion orientation warnings without running conversion."""

import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from conv2Nifti_auto import warn_if_niftis_not_lip


class ConversionOrientationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def save(self, relative, orientation, shape=(2, 3, 4)):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        affine = np.eye(4)
        affine[:3, :] = 0
        vectors = {"L": (-1, 0, 0), "R": (1, 0, 0), "I": (0, 0, -1),
                   "S": (0, 0, 1), "P": (0, -1, 0), "A": (0, 1, 0)}
        for axis, code in enumerate(orientation):
            affine[:3, axis] = vectors[code]
        nib.save(nib.Nifti1Image(np.zeros(shape), affine), path)
        return path

    def check_output(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            warn_if_niftis_not_lip(self.root)
        return output.getvalue()

    def test_lip_3d_and_4d_outputs_and_empty_folder_are_silent(self):
        self.assertEqual(self.check_output(), "")
        self.save("sub-01/ses-01/anat/lip.nii", "LIP")
        self.save("sub-01/ses-01/dwi/derived/lip.nii.gz", "LIP", (2, 3, 4, 2))
        self.assertEqual(self.check_output(), "")

    def test_non_lip_files_warn_with_helper_path_and_remain_unchanged(self):
        self.save("sub-01/ses-01/anat/lip.nii", "LIP")
        ras = self.save("sub-01/ses-01/func/ras.nii.gz", "RAS", (2, 3, 4, 2))
        rpi = self.save("sub-01/ses-01/t2map/rpi.nii", "RPI")
        before = {p: p.read_bytes() for p in (ras, rpi)}
        with self.assertLogs(level="WARNING") as logs:
            output = self.check_output()
        self.assertIn("2 NIfTI file(s) are not in LIP orientation", output)
        self.assertIn("bin/helper_tools/ReorientBatch.py", output)
        self.assertIn("ras.nii.gz: RAS", output)
        self.assertIn("rpi.nii: RPI", output)
        self.assertNotIn("lip.nii:", output)
        self.assertTrue(any("RAS" in line for line in logs.output))
        for path, contents in before.items():
            self.assertEqual(path.read_bytes(), contents)

    def test_fieldmap_and_magnitude_names_are_excluded_from_orientation_check(self):
        for name in ("sub-01_fieldmap.nii", "sub-01_magnitude.nii.gz",
                     "sub-01_FieldMap_extra.nii.gz", "sub-01_MAGNITUDE_extra.nii"):
            self.save(f"sub-01/ses-01/fmap/{name}", "IPL")
        with self.assertNoLogs(level="WARNING"):
            self.assertEqual(self.check_output(), "")
        # A folder name alone must not exclude an unrelated NIfTI.
        self.save("fieldmap/sub-01_anat.nii.gz", "RAS")
        with self.assertLogs(level="WARNING"):
            output = self.check_output()
        self.assertIn("1 NIfTI file(s) are not in LIP orientation", output)
        self.assertIn("sub-01_anat.nii.gz: RAS", output)

    def test_unreadable_header_does_not_hide_orientation_warning(self):
        (self.root / "broken.nii.gz").write_bytes(b"not a nifti")
        (self.root / "ignore.json").write_text("{}")
        self.save("ras.nii", "RAS")
        with self.assertLogs(level="WARNING"):
            output = self.check_output()
        self.assertIn("1 NIfTI file(s) are not in LIP orientation", output)
        self.assertIn("Orientation could not be checked for 1 NIfTI file(s)", output)
        self.assertIn("broken.nii.gz", output)


if __name__ == "__main__":
    unittest.main()
