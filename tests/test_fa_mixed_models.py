"""Regression checks for both standalone FA analyses and their primary models."""

import csv
import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load_analysis(region):
    path = ROOT / f"bin/3.2_DTIConnectivity/{region}_fa_group_analysis.py"
    spec = importlib.util.spec_from_file_location(f"{region}_mixed_tests", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_data():
    rng = np.random.default_rng(721)
    rows = []
    for animal in range(32):
        group = "Sham" if animal < 16 else "tDCS"
        intercept = rng.normal(0, 0.04)
        for timepoint in ("PT3", "PT14"):
            lesion = rng.uniform(2, 18)
            for side in ("left", "right"):
                value = (0.45 + intercept + 0.015 * (group == "tDCS")
                         + 0.02 * (timepoint == "PT14") - 0.08 * (side == "left")
                         + 0.002 * lesion + rng.normal(0, 0.009))
                rows.append(dict(subject=f"sub-{animal:03}", group=group,
                                 timepoint=timepoint, side=side, median_fa=value,
                                 lesion_volume_mm3=lesion))
    data = pd.DataFrame(rows)
    # Whole missing sessions and incomplete sides remain in the long QC table.
    data.loc[(data.subject == "sub-000") & (data.timepoint == "PT14"), "median_fa"] = np.nan
    data.loc[(data.subject == "sub-001") & (data.timepoint == "PT3"), "median_fa"] = np.nan
    data.loc[9, "median_fa"] = np.inf
    data.loc[data.subject == "sub-003", "lesion_volume_mm3"] = np.nan
    data.loc[(data.subject == "sub-004") & (data.timepoint == "PT3"), "lesion_volume_mm3"] = np.inf
    return data


class MixedAnalysisChecks:
    def setUp(self):
        self.analysis = load_analysis(self.region)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def save(self, path, data, affine=None, units="mm"):
        path.parent.mkdir(parents=True, exist_ok=True)
        image = nib.Nifti1Image(np.asarray(data, dtype=float), np.eye(4) if affine is None else affine)
        image.header.set_xyzt_units(units)
        nib.save(image, path)
        return path

    def test_physical_volume_and_units(self):
        mask = np.zeros((2, 2, 2))
        mask.flat[:3] = 1
        for units, factor in (("mm", 1), ("meter", .001), ("micron", 1000)):
            with self.subTest(units=units):
                # Reflection plus shear: determinant, not a voxel-count surrogate.
                affine = np.array([[-.2, .1, 0, 0], [0, .3, 0, 0], [0, 0, .4, 0], [0, 0, 0, 1.]])
                affine[:3, :3] *= factor
                path = self.save(self.folder / f"mask_{units}.nii.gz", mask, affine, units)
                volume = self.analysis.calculate_lesion_volume(path)
                self.assertEqual(volume["n_lesion_voxels"], 3)
                self.assertAlmostEqual(volume["lesion_volume_mm3"], .072, places=7)
        empty = self.save(self.folder / "empty.nii", np.zeros((2, 2, 2)))
        self.assertEqual(self.analysis.calculate_lesion_volume(empty)["lesion_volume_mm3"], 0)

    def test_invalid_mask_qc(self):
        for invalid in (np.nan, np.inf, .5, 2, -1):
            with self.subTest(value=invalid):
                mask = np.zeros((2, 2, 2))
                mask.flat[0] = invalid
                path = self.save(self.folder / "bad.nii", mask)
                with self.assertRaisesRegex(ValueError, "finite numeric 0 and 1"):
                    self.analysis.calculate_lesion_volume(path)
        for shape in ((2, 2), (2, 2, 2, 1)):
            path = self.save(self.folder / "dimensions.nii", np.zeros(shape))
            with self.assertRaisesRegex(ValueError, "3D"):
                self.analysis.calculate_lesion_volume(path)
        path = self.save(self.folder / "unknown.nii", np.ones((2, 2, 2)), units="unknown")
        with self.assertRaisesRegex(ValueError, "spatial units"):
            self.analysis.calculate_lesion_volume(path)
        for affine in (np.diag([0, 1, 1, 1]), np.full((4, 4), np.nan)):
            with patch.object(self.analysis.nib, "load", return_value=SimpleNamespace(shape=(2, 2, 2), affine=affine)):
                with self.assertRaises(ValueError):
                    self.analysis.calculate_lesion_volume(path)

    def test_session_masks_long_format_and_ambiguity(self):
        rows = []
        for subject in ("sub-A", "sub-B"):
            for timepoint in self.analysis.TIMEPOINTS:
                for side in self.analysis.LABELS:
                    rows.append(dict(subject=subject, group="Sham", timepoint=timepoint,
                                     side=side, median_fa=.4, fa_file="fa.nii", atlas_file="atlas.nii", status="ok"))
        rows[1]["median_fa"] = np.nan
        for timepoint, count in (("PT3", 2), ("PT14", 5)):
            mask = np.zeros((2, 2, 2))
            mask.flat[:count] = 1
            self.save(self.folder / f"sub-A/ses-{timepoint}/dwi/scanBetStroke_mask.nii.gz", mask)
        ambiguous = self.folder / "sub-B/ses-PT3/dwi"
        for suffix in ("nii", "nii.gz"):
            self.save(ambiguous / f"scanBetStroke_mask.{suffix}", np.ones((2, 2, 2)))
        with self.assertRaisesRegex(ValueError, "Ambiguous lesion"):
            self.analysis.find_lesion_mask(ambiguous)
        with self.assertLogs(self.analysis.LOG, level="INFO") as logs:
            data = self.analysis.build_long_format(rows, self.folder)
        self.assertIn("FA present=True, atlas present=True", "\n".join(logs.output))
        self.assertEqual(data.included_model_1.sum(), 7)
        self.assertEqual(data.included_model_2.sum(), 3)
        self.assertEqual(list(data.lesion_volume_mm3[:4]), [2, 2, 5, 5])
        self.assertEqual(data.lesion_status[4], "ambiguous_mask")
        self.assertEqual(data.lesion_status[6], "missing_mask")
        np.testing.assert_allclose(data.lesion_volume_centered[[0, 2, 3]], [-2, 1, 1])
        self.assertEqual(rows[0]["side"], "ipsilesional")
        self.assertEqual(data.side[0], "left")
        # Invalid masks stay local to model 2 and retain the chosen path.
        path = self.save(self.folder / "sub-B/ses-PT14/dwi/scanBetStroke_mask.nii", np.full((2, 2, 2), .4))
        data = self.analysis.build_long_format(rows, self.folder)
        self.assertEqual(data.lesion_status[6], "invalid_mask")
        self.assertEqual(data.lesion_mask_file[6], str(path))
        self.assertTrue(data.included_model_1[6])
        self.assertFalse(data.included_model_2[6])

    def test_real_unbalanced_models_and_exports(self):
        data = synthetic_data()
        original = data.copy(deep=True)
        for with_lesion in (False, True):
            with self.subTest(with_lesion=with_lesion):
                coefficients, meta = self.analysis.fit_mixed_model(data, with_lesion)
                expected = np.isfinite(data.median_fa)
                if with_lesion:
                    expected &= np.isfinite(data.lesion_volume_mm3)
                self.assertTrue(meta["converged"], meta)
                self.assertEqual(meta["status"], "ok", meta)
                self.assertEqual(meta["n_subjects"], data.loc[expected, "subject"].nunique())
                self.assertEqual(meta["n_observations"], expected.sum())
                self.assertEqual(len(coefficients), 9 if with_lesion else 8)
                self.assertEqual(meta["reference_side"], "right")
                indexed = {c["term"]: c for c in coefficients}
                self.assertAlmostEqual(indexed["C(side)[T.left]"]["estimate"], -.08, delta=.02)
                self.assertTrue(np.isfinite([meta[k] for k in ("AIC", "BIC", "log_likelihood", "random_effect_variance")]).all())
                if with_lesion:
                    self.assertAlmostEqual(meta["lesion_centering_mean_mm3"], data.loc[expected, "lesion_volume_mm3"].mean())
                    self.assertAlmostEqual(indexed["lesion_volume_centered"]["estimate"], .002, delta=.001)
                stem = "test_with_lesion" if with_lesion else "test"
                self.analysis.write_mixed_model_results(self.folder, stem, coefficients, meta)
                self.assertTrue((self.folder / f"{stem}.csv").is_file())
                self.assertIn("reference_side: right", (self.folder / f"{stem}_summary.txt").read_text())
        pd.testing.assert_frame_equal(data, original)

        # Even an animal with just one finite observation stays in model 1.
        singleton = data.index[data.subject == "sub-005"]
        data.loc[singleton[:3], "median_fa"] = np.nan
        _, meta = self.analysis.fit_mixed_model(data)
        self.assertEqual(meta["n_subjects"], 32)
        self.assertEqual(meta["n_observations"], 120)
        self.assertTrue(meta["converged"])

    def test_unidentifiable_and_failed_models_report_all_terms(self):
        data = synthetic_data()
        cases = []
        no_lesion = data.copy()
        no_lesion["lesion_volume_mm3"] = np.nan
        cases.append((no_lesion, True, "insufficient_data"))
        missing_cell = data.loc[~((data.group == "tDCS") & (data.side == "left") & (data.timepoint == "PT14"))]
        cases.append((missing_cell, False, "rank_deficient"))
        constant_lesion = data.copy()
        constant_lesion["lesion_volume_mm3"] = 3
        cases.append((constant_lesion, True, "rank_deficient"))
        singletons = data.copy()
        singletons["subject"] = [f"sub-{i}" for i in range(len(data))]
        cases.append((singletons, False, "no_repeated_observations"))
        for frame, with_lesion, status in cases:
            coefficients, meta = self.analysis.fit_mixed_model(frame, with_lesion)
            self.assertEqual(meta["status"], status)
            self.assertFalse(meta["converged"])
            self.assertEqual(len(coefficients), 9 if with_lesion else 8)
            self.assertTrue(all(np.isnan(c["estimate"]) for c in coefficients))
        with patch.object(self.analysis.smf, "mixedlm") as constructor:
            constructor.return_value.exog = np.eye(8)
            constructor.return_value.fit.side_effect = np.linalg.LinAlgError("synthetic failure")
            coefficients, meta = self.analysis.fit_mixed_model(data)
            self.assertEqual(meta["status"], "fit_failed")
            self.assertIn("synthetic failure", meta["warnings"])
            self.assertEqual(constructor.return_value.fit.call_count, 2)

    def test_nonconvergence_is_exported_and_warned(self):
        data = synthetic_data()
        finite = data.loc[np.isfinite(data.median_fa)].copy()
        for col, categories in self.analysis.MODEL_CATEGORIES.items():
            finite[col] = pd.Categorical(finite[col], categories=categories)
        model = self.analysis.smf.mixedlm(self.analysis.MIXED_FORMULA, finite,
                                        groups=finite.subject, re_formula="1")
        result = model.fit(reml=False, method="powell", disp=False)
        result.converged = False
        with patch.object(self.analysis.smf, "mixedlm", return_value=model), patch.object(model, "fit", return_value=result):
            with self.assertLogs(self.analysis.LOG, level="WARNING"):
                coefficients, meta = self.analysis.fit_mixed_model(data)
        self.assertEqual(meta["status"], "not_converged")
        self.assertFalse(meta["converged"])
        self.assertTrue(np.isfinite(coefficients[0]["estimate"]))
        self.analysis.write_mixed_model_results(self.folder, "nonconverged", coefficients, meta)
        self.assertIn("converged: False", (self.folder / "nonconverged_summary.txt").read_text())

    def test_end_to_end_retains_old_outputs_and_adds_new_outputs(self):
        data = synthetic_data()
        groups = self.folder / "groups.csv"
        groups.write_text("tDCS,sham\n" + "".join(f"sub-{i+16:03},sub-{i:03}\n" for i in range(16)))
        for (subject, timepoint), session in data.groupby(["subject", "timepoint"]):
            if not np.isfinite(session.median_fa).any():
                continue
            dwi = self.folder / subject / f"ses-{timepoint}" / "dwi"
            labels = np.zeros((2, 2, 2))
            labels[0], labels[1] = self.analysis.LABELS.values()
            values = session.set_index("side").median_fa
            fa = np.empty((2, 2, 2))
            fa[0], fa[1] = values["left"], values["right"]
            self.save(dwi / "DSI_studio/scan.fa.nii.gz", fa)
            atlas_name = "scanAnnoSplit.nii.gz" if self.region == "ic" else "scan_AnnoSplit_parental.nii.gz"
            self.save(dwi / atlas_name, labels)
            lesion_volume = session.lesion_volume_mm3.iloc[0]
            if np.isfinite(lesion_volume):
                affine = np.diag([lesion_volume / 8, 1, 1, 1])
                self.save(dwi / "scanBetStroke_mask.nii.gz", np.ones((2, 2, 2)), affine)
        output = self.folder / "output"
        self.assertEqual(self.analysis.main(["-i", str(self.folder), "-g", str(groups), "-o", str(output)]), 0)
        prefix = f"{self.region}_fa"
        for name in ("animal_medians.csv", "statistics.csv", "stroke_effect_statistics.csv", "long_format.csv",
                     "mixed_model.csv", "mixed_model_with_lesion.csv", "mixed_model_summary.txt", "mixed_model_with_lesion_summary.txt"):
            self.assertGreater((output / f"{prefix}_{name}").stat().st_size, 100)
        for name in ("PT3", "PT14", "stroke_effect_PT3", "stroke_effect_PT14", "all_values_mixed_model"):
            for extension in ("png", "pdf"):
                self.assertGreater((output / f"{prefix}_{name}.{extension}").stat().st_size, 1000)
        long_data = pd.read_csv(output / f"{prefix}_long_format.csv")
        self.assertEqual(len(long_data), 128)
        self.assertEqual(long_data.included_model_1.sum(), 123)
        self.assertEqual(long_data.included_model_2.sum(), 117)
        for lesion_suffix, expected_n in (("", 123), ("_with_lesion", 117)):
            stats = pd.read_csv(output / f"{prefix}_mixed_model{lesion_suffix}.csv")
            self.assertTrue((stats.n_observations == expected_n).all())
        # A completely empty FA dataset still writes model failure reports.
        empty_output = self.folder / "empty_output"
        empty_proc = self.folder / "empty_proc"
        empty_proc.mkdir()
        self.assertEqual(self.analysis.main(["-i", str(empty_proc), "-g", str(groups), "-o", str(empty_output)]), 1)
        with (empty_output / f"{prefix}_mixed_model.csv").open() as handle:
            self.assertEqual(next(csv.DictReader(handle))["status"], "insufficient_data")


class ICMixedTests(MixedAnalysisChecks, unittest.TestCase):
    region = "ic"


class CCMixedTests(MixedAnalysisChecks, unittest.TestCase):
    region = "cc"


if __name__ == "__main__":
    unittest.main()
