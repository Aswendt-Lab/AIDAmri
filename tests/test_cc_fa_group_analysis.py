"""Check CC extraction, spatial alignment and animal-level comparisons."""

import importlib.util
from pathlib import Path
import tempfile
import unittest

import nibabel as nib
import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "bin/3.2_DTIConnectivity/cc_fa_group_analysis.py"
spec = importlib.util.spec_from_file_location("cc_fa_group_analysis", SCRIPT)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class CCAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def save(self, name, data, affine=None):
        path = self.folder / name
        path.parent.mkdir(parents=True, exist_ok=True)
        nib.save(nib.Nifti1Image(np.asarray(data, dtype=float),
                               np.eye(4) if affine is None else affine), path)
        return path

    def test_wide_csv_and_conflicting_assignment(self):
        path = self.folder / "groups.csv"
        path.write_text("tDCS,sham\nsub-A,sub-B\n,sub-C\n")
        self.assertEqual(analysis.load_groups(path), {"sub-A": "tDCS", "sub-B": "Sham", "sub-C": "Sham"})
        path.write_text("tDCS,sham\nsub-A,sub-A\n")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            analysis.load_groups(path)

    def test_zero_nan_and_infinity_excluded_before_median(self):
        fa = self.save("fa.nii.gz", np.array([0, .2, .4, np.nan, .5, .9, np.inf, 0]).reshape(2, 2, 2))
        atlas = self.save("atlas.nii.gz", np.array([776]*4 + [2776]*4).reshape(2, 2, 2))
        left, right = analysis.extract_cc(fa, atlas)
        self.assertAlmostEqual(left["median_fa"], .3)
        self.assertAlmostEqual(right["median_fa"], .7)
        for row in (left, right):
            self.assertEqual(row["n_valid_voxels"], 2)
            self.assertEqual(row["n_zero_voxels"], 1)
            self.assertEqual(row["n_nonfinite_voxels"], 1)

    def test_affine_flip_preserves_left_right_anatomy(self):
        values = np.zeros((4, 2, 2))
        values[0] = .2
        values[3] = .8
        labels = np.zeros_like(values)
        labels[3] = 776
        labels[0] = 2776
        affine = np.eye(4)
        affine[0, 0], affine[0, 3] = -1, 3
        fa = self.save("fa.nii", values)
        atlas = self.save("atlas.nii", labels, affine)
        left, right = analysis.extract_cc(fa, atlas)
        self.assertTrue(left["atlas_resampled"])
        self.assertAlmostEqual(left["median_fa"], .2)
        self.assertAlmostEqual(right["median_fa"], .8)

    def test_empty_and_all_zero_rois_remain_nan(self):
        fa = self.save("fa.nii", np.zeros((2, 2, 2)))
        atlas = self.save("atlas.nii", np.full((2, 2, 2), 776))
        left, right = analysis.extract_cc(fa, atlas)
        self.assertEqual(left["status"], "no_valid_fa")
        self.assertEqual(right["status"], "empty_roi")
        self.assertTrue(np.isnan(left["median_fa"]))
        self.assertTrue(np.isnan(right["median_fa"]))

    def test_holm_keeps_four_planned_comparisons(self):
        np.testing.assert_allclose(analysis.holm_adjust([.01, .04, .03, .2]), [.04, .09, .09, .2])
        np.testing.assert_allclose(analysis.holm_adjust([.01, np.nan, .04, np.nan]),
                                   [.04, np.nan, .12, np.nan], equal_nan=True)

    def test_missing_fa_reported_and_ambiguous_fa_rejected(self):
        dwi = self.folder / "sub-A/ses-PT3/dwi"
        dwi.mkdir(parents=True)
        rows = analysis.collect_medians(self.folder, {"sub-A": "tDCS"})
        self.assertEqual(rows[0]["status"], "missing_fa")
        self.assertEqual(rows[2]["status"], "missing_session")
        for name in ("a.fa.nii", "b.fa0.nii.gz"):
            self.save(f"sub-A/ses-PT3/dwi/DSI_studio/{name}", np.zeros((2, 2, 2)))
        with self.assertRaisesRegex(ValueError, "Ambiguous FA"):
            analysis.find_inputs(dwi)

    def test_end_to_end_has_one_observation_per_animal_side_timepoint(self):
        groups = self.folder / "groups.csv"
        groups.write_text("tDCS,sham\nsub-A,sub-C\nsub-B,sub-D\n")
        labels = np.zeros((2, 2, 2))
        labels[0], labels[1] = 776, 2776
        for index, subject in enumerate(("A", "B", "C", "D")):
            for tp in analysis.TIMEPOINTS:
                dwi = f"sub-{subject}/ses-{tp}/dwi"
                self.save(f"{dwi}/DSI_studio/scan.fa.nii.gz", np.full(labels.shape, (index+1)/10))
                self.save(f"{dwi}/scan_AnnoSplit_parental.nii.gz", labels)
        output = self.folder / "results"
        self.assertEqual(analysis.main(["-i", str(self.folder), "-g", str(groups), "-o", str(output)]), 0)
        import csv
        with (output / "cc_fa_animal_medians.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 16)
        with (output / "cc_fa_statistics.csv").open() as handle:
            stats = list(csv.DictReader(handle))
        self.assertEqual(len(stats), 4)
        self.assertTrue(all(s["n_sham"] == s["n_tdcs"] == "2" for s in stats))
        for tp in analysis.TIMEPOINTS:
            for extension in ("png", "pdf"):
                self.assertGreater((output / f"cc_fa_{tp}.{extension}").stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
