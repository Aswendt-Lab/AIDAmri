"""Check MATLAB exports preserve DSI data and connectivity definitions."""

import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.io import loadmat, savemat

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin/3.2_DTIConnectivity"))
import dsi_tools


def fixture(value=7):
    counts = np.array([[value], [5]], dtype=np.float32)
    return {
        "connectivity": np.array([[0, 3], [3, 0]], dtype=np.float32),
        "number of tracts t2r": counts,
        "number of tracts r2r": np.array([[0, 3], [3, 0]], dtype=np.float32),
        "mean length(mm) t2r": np.array([[1.25], [2.5]], dtype=np.float32),
        "name": np.frombuffer(b"L_Primary_motor_area\nR_Primary_motor_area\n\0", dtype=np.uint8)[None, :],
        "atlas": np.frombuffer(b"test_atlas\0", dtype=np.uint8)[None, :],
    }


class MatlabExportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.path = self.folder / "result.mat"

    def test_v4_metrics_labels_and_legacy_fields_survive_conversion(self):
        original = fixture()
        savemat(self.path, original, format="4")
        self.path.chmod(0o640)
        dsi_tools.make_connectivity_matlab_compatible(self.path)
        result = loadmat(self.path)
        for key in ("connectivity", "name", "atlas"):
            np.testing.assert_array_equal(result[key], original[key])
        np.testing.assert_array_equal(result["number_of_tracts_t2r"], original["number of tracts t2r"])
        np.testing.assert_array_equal(result["number_of_tracts_r2r"], original["number of tracts r2r"])
        np.testing.assert_array_equal(result["mean_length_mm_t2r"], original["mean length(mm) t2r"])
        self.assertEqual(result["number_of_tracts_t2r"].dtype, np.float32)
        self.assertEqual([x[0] for x in result["region_names"].ravel()],
                         ["L_Primary_motor_area", "R_Primary_motor_area"])
        pairs = [(left[0], right[0]) for left, right in result["matlab_variable_name_map"]]
        self.assertIn(("mean length(mm) t2r", "mean_length_mm_t2r"), pairs)
        self.assertTrue(self.path.read_bytes().startswith(b"MATLAB 5.0 MAT-file"))
        self.assertEqual(self.path.stat().st_mode & 0o777, 0o640)

    def test_name_collisions_and_long_identifiers_do_not_lose_metrics(self):
        original = {"fa t2r": np.array([[1]]), "fa-t2r": np.array([[2]]),
                    "fa_t2r": np.array([[3]]), "9 tracts": np.array([[4]]),
                    "a" * 64: np.array([[5]]), "a" * 63 + "b": np.array([[6]])}
        savemat(self.path, original, format="4")
        dsi_tools.make_connectivity_matlab_compatible(self.path)
        result = loadmat(self.path)
        names = []
        for left, right in result["matlab_variable_name_map"]:
            old, new = left[0], right[0]
            names.append(new)
            self.assertRegex(new, r"^[A-Za-z][A-Za-z0-9_]{0,62}$")
            np.testing.assert_array_equal(result[new], original[old])
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual(int(result["fa_t2r"][0, 0]), 3)

    def test_repeated_conversion_preserves_existing_mapping(self):
        savemat(self.path, fixture(), format="4")
        dsi_tools.make_connectivity_matlab_compatible(self.path)
        first = loadmat(self.path)
        dsi_tools.make_connectivity_matlab_compatible(self.path)
        second = loadmat(self.path)
        self.assertEqual(set(first), set(second))
        for key in ("number_of_tracts_t2r", "region_names", "matlab_variable_name_map"):
            np.testing.assert_equal(first[key], second[key])

    def test_failed_save_preserves_original_and_removes_temporary_file(self):
        savemat(self.path, fixture(), format="4")
        original = self.path.read_bytes()
        with patch.object(dsi_tools, "savemat", side_effect=OSError("Disk full")), \
                self.assertRaisesRegex(OSError, "Disk full"):
            dsi_tools.make_connectivity_matlab_compatible(self.path)
        self.assertEqual(self.path.read_bytes(), original)
        self.assertEqual(list(self.folder.iterdir()), [self.path])

    def test_pipeline_converts_pass_and_end_exports_before_moving(self):
        fib_dir = self.folder / "fib_map"
        fib_dir.mkdir()
        (fib_dir / "scan.fz").touch()
        (fib_dir / "scan.fz.trk.gz").touch()
        seed = self.folder / "atlas.nii.gz"
        seed.touch()

        def run(command, **kwargs):
            kind = next(x.split("=", 1)[1] for x in command if x.startswith("--connectivity_type="))
            savemat(fib_dir / "scan.fz.trk.gz.atlas.connectivity.mat",
                    fixture(7 if kind == "pass" else 4), format="4")

        with patch.object(dsi_tools.subprocess, "run", side_effect=run), \
                patch.object(dsi_tools, "create_stroke_mask_label_file"), \
                contextlib.redirect_stdout(io.StringIO()):
            dsi_tools.connectivity("dsi_studio", str(fib_dir), str(seed),
                                   str(self.folder), "connectivity")
        for kind, expected in (("pass", 7), ("end", 4)):
            path = self.folder / "connectivity" / f"scan.fz.trk.gz.atlas.all.{kind}.connectivity.mat"
            result = loadmat(path)
            self.assertEqual(int(result["number_of_tracts_t2r"][0, 0]), expected)
            self.assertIn("region_names", result)
        self.assertEqual(list(fib_dir.glob("*.mat")), [])
        self.assertEqual(list((self.folder / 'connectivity').glob('*.csv')), [])


if __name__ == "__main__":
    unittest.main()
