"""Validate subject-level models, contrasts, correction and plot annotations."""

import contextlib
import csv
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np
from scipy import stats as distributions
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[1]


def load_script(name, filename):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin/3.1_T2Processing" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


analysis = load_script("affected_regions_statistics", "affected_regions_statistics.py")
plots = load_script("affected_regions_plots", "affected_regions_plots.py")
aggregate = load_script("statistical_cohort_export", "aggregate_affected_regions.py")
OUTCOME = "L_Region [%]"


def synthetic_rows(boundary=False):
    rng = np.random.default_rng(142)
    rows = []
    count = 12 if boundary else 20
    for group in ("Sham", "tDCS"):
        for subject in range(count):
            animal = 0 if boundary else rng.normal(0, 1.5)
            for time, session in enumerate(("ses-PT3", "ses-PT14")):
                error = ((subject - 5.5) / 6 * (1 if time == 0 else -1)
                         if boundary else rng.normal(0, .4))
                value = 10 + animal + 2 * (group == "tDCS") + time + 3 * (group == "tDCS") * time + error
                rows.append(dict(subject=f"{group}-{subject}", group=group, session=session,
                                 **{OUTCOME: value, "L_Region [mm3]": value / 5,
                                    "Total stroke volume [mm3]": value * 2,
                                    "Stroke percent of brain [%]": value / 2}))
    return rows


def cell_design(rows):
    cells = [(group, session) for group in ("Sham", "tDCS") for session in ("ses-PT3", "ses-PT14")]
    matrix = np.array([[float((row["group"], row["session"]) == cell) for cell in cells] for row in rows])
    subjects = sorted({row["subject"] for row in rows})
    z = np.array([[float(row["subject"] == subject) for subject in subjects] for row in rows])
    return cells, matrix, z


class AffectedRegionStatisticsTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)

    def test_reml_contrasts_match_independent_generalized_least_squares(self):
        rows = synthetic_rows()
        result = analysis.analyze_endpoint(rows, OUTCOME)
        meta = result["metadata"]
        self.assertEqual(meta["status"], "ok", meta)
        self.assertEqual((meta["n_subjects"], meta["n_observations"], meta["family_size"]), (40, 80, 4))
        self.assertEqual(meta["fit_method"], "REML MixedLM")
        cells, x, z = cell_design(rows)
        variance = meta["random_intercept_variance"] * z @ z.T + meta["residual_variance"] * np.eye(len(rows))
        inverse_x = np.linalg.solve(variance, x)
        covariance = np.linalg.inv(x.T @ inverse_x)
        beta = covariance @ x.T @ np.linalg.solve(variance, np.array([row[OUTCOME] for row in rows]))
        for contrast in result["contrasts"]:
            difference = np.zeros(4)
            difference[cells.index((contrast["group_a"], contrast["session_a"]))] = 1
            difference[cells.index((contrast["group_b"], contrast["session_b"]))] = -1
            estimate = float(difference @ beta)
            se = float(np.sqrt(difference @ covariance @ difference))
            self.assertAlmostEqual(contrast["estimate"], estimate, places=6)
            self.assertAlmostEqual(contrast["standard_error"], se, places=5)
            self.assertAlmostEqual(contrast["p_raw"], 2 * distributions.norm.sf(abs(estimate / se)), places=6)
            self.assertEqual(contrast["test"], "z")
        effects = {row["term"]: row for row in result["tests"]}
        self.assertEqual(set(effects), {"group", "time", "interaction"})
        self.assertLess(effects["interaction"]["p_raw"], .001)
        self.assertTrue(all(row["test"] == "chi-square" for row in effects.values()))

    def test_holm_sidak_is_monotone_and_corrected_separately_per_outcome(self):
        rows = synthetic_rows()
        with contextlib.redirect_stdout(io.StringIO()):
            results = analysis.analyze_statistics(rows, [OUTCOME, "Total stroke volume [mm3]"])
        for result in results.values():
            contrasts = result["contrasts"]
            self.assertEqual(len(contrasts), 4)
            ordered = sorted(contrasts, key=lambda row: row["p_raw"])
            previous = 0
            for rank, row in enumerate(ordered):
                corrected = max(previous, -np.expm1((4 - rank) * np.log1p(-row["p_raw"])))
                self.assertAlmostEqual(row["p_adjusted"], corrected)
                self.assertEqual(row["family_size"], 4)
                previous = corrected
        np.testing.assert_allclose([row["p_adjusted"] for row in results[OUTCOME]["contrasts"]],
                                   [row["p_adjusted"] for row in results["Total stroke volume [mm3]"]["contrasts"]],
                                   rtol=.001)

    def test_more_than_two_groups_and_timepoints_use_all_simple_contrasts(self):
        rng = np.random.default_rng(91)
        rows = []
        for group in range(3):
            for subject in range(10):
                offset = rng.normal(0, 1)
                for time in range(3):
                    rows.append(dict(subject=f"{group}-{subject}", group=f"group-{group}", session=f"ses-PT{time}",
                                     **{OUTCOME: 10 + group * 2 + time + group * time + offset + rng.normal(0, .3)}))
        result = analysis.analyze_endpoint(rows, OUTCOME)
        self.assertEqual(result["metadata"]["status"], "ok", result["metadata"])
        self.assertEqual(len(result["contrasts"]), 18)
        self.assertTrue(all(row["family_size"] == 18 for row in result["contrasts"]))
        self.assertEqual([row["df_num"] for row in result["tests"]], [2, 2, 4])

    def test_zero_random_variance_refits_full_ols_and_matches_classical_t_tests(self):
        rows = synthetic_rows(boundary=True)
        result = analysis.analyze_endpoint(rows, OUTCOME)
        meta = result["metadata"]
        self.assertEqual(meta["status"], "simplified_ols", meta)
        self.assertLessEqual(meta["random_variance_ratio"], 1e-8)
        self.assertEqual(meta["random_intercept_variance"], 0)
        cells, x, _ = cell_design(rows)
        reference = sm.OLS([row[OUTCOME] for row in rows], x).fit()
        self.assertAlmostEqual(meta["residual_variance"], reference.scale)
        for row in result["contrasts"]:
            c = np.zeros(4)
            c[cells.index((row["group_a"], row["session_a"]))] = 1
            c[cells.index((row["group_b"], row["session_b"]))] = -1
            test = reference.t_test(c)
            self.assertAlmostEqual(row["estimate"], float(np.asarray(test.effect).item()))
            self.assertAlmostEqual(row["standard_error"], float(np.asarray(test.sd).item()))
            self.assertAlmostEqual(row["p_raw"], float(np.asarray(test.pvalue).item()))
            self.assertEqual(row["df_denom"], reference.df_resid)
        self.assertTrue(all(row["test"] == "F" for row in result["tests"]))

    def test_measurement_rescaling_does_not_change_tests_or_boundary_decision(self):
        for boundary in (False, True):
            with self.subTest(boundary=boundary):
                rows = synthetic_rows(boundary=boundary)
                first = analysis.analyze_endpoint(rows, OUTCOME)
                scaled = [dict(row, **{OUTCOME: row[OUTCOME] * 1e-6}) for row in rows]
                second = analysis.analyze_endpoint(scaled, OUTCOME)
                self.assertEqual(first["metadata"]["status"], second["metadata"]["status"])
                np.testing.assert_allclose([row["p_adjusted"] for row in first["contrasts"]],
                                           [row["p_adjusted"] for row in second["contrasts"]], rtol=1e-4)
                for a, b in zip(first["contrasts"], second["contrasts"]):
                    self.assertAlmostEqual(a["estimate"] * 1e-6, b["estimate"], places=10)

    def test_zeros_remain_missing_values_drop_and_unassigned_subjects_are_excluded(self):
        rows = synthetic_rows()
        rows[0][OUTCOME] = 0
        rows[1][OUTCOME] = None
        rows[2][OUTCOME] = np.nan
        rows.extend(dict(row, subject="unknown", group="Unassigned") for row in rows[:2])
        original = [dict(row) for row in rows]
        result = analysis.analyze_endpoint(rows, OUTCOME)
        meta = result["metadata"]
        self.assertIn(meta["status"], ("ok", "simplified_ols"), meta)
        self.assertEqual(meta["n_observations"], 78)
        self.assertEqual(meta["n_zero_values"], 1)
        self.assertEqual(meta["n_missing_values"], 2)
        self.assertEqual(meta["n_unassigned_observations"], 2)
        self.assertEqual(rows, original)

    def test_invalid_or_unidentifiable_models_have_no_pvalues_or_stars(self):
        rows = synthetic_rows()
        cases = []
        small = [row for row in rows if row["subject"] in ("Sham-0", "tDCS-0")]
        cases.append((small, "insufficient_data"))
        cases.append(([dict(row, **{OUTCOME: 0}) for row in rows], "constant_data"))
        cases.append(([row for row in rows if not (row["group"] == "tDCS" and row["session"] == "ses-PT14")],
                      "rank_deficient"))
        cases.append(([dict(row, subject=f"single-{index}") for index, row in enumerate(rows)], "insufficient_data"))
        changed = [dict(row) for row in rows]
        changed[1]["group"] = "tDCS"
        cases.append((changed, "invalid_data"))
        cases.append(([dict(row, **{OUTCOME: None}) for row in rows], "insufficient_data"))
        for data, expected_status in cases:
            with self.subTest(status=expected_status):
                result = analysis.analyze_endpoint(data, OUTCOME)
                self.assertEqual(result["metadata"]["status"], expected_status, result["metadata"])
                self.assertTrue(all(row["p_raw"] is None for row in result["tests"]))
                self.assertTrue(all(row["p_raw"] is None and row["p_adjusted"] is None and not row["stars"]
                                    for row in result["contrasts"]))

    def test_nonconvergence_or_invalid_covariance_does_not_trigger_ols(self):
        invalid = SimpleNamespace(converged=True, cov_re=np.array([[1.]]), scale=1., llf=-10,
                                  fe_params=np.zeros(4), cov_params=lambda: np.diag([1., -1., 1., 1.]))
        for fitted, expected in ((SimpleNamespace(converged=False), "not_converged"), (invalid, "fit_failed")):
            with self.subTest(expected=expected), mock.patch.object(analysis, "MixedLM") as constructor, \
                    mock.patch.object(analysis.sm, "OLS") as ols:
                constructor.return_value.fit.return_value = fitted
                result = analysis.analyze_endpoint(synthetic_rows(), OUTCOME)
                self.assertEqual(result["metadata"]["status"], expected)
                self.assertTrue(all(row["p_adjusted"] is None for row in result["contrasts"]))
                ols.assert_not_called()

    def test_significance_thresholds_use_strict_adjusted_pvalues(self):
        for p, expected in ((None, ""), (np.nan, ""), (.05, ""), (.049, "*"), (.01, "*"),
                            (.009, "**"), (.001, "**"), (.0009, "***"), (.0001, "***"), (.00009, "****")):
            self.assertEqual(analysis.significance_stars(p), expected)

    def test_plot_annotations_use_only_adjusted_significance_and_do_not_label_failed_tests_ns(self):
        rows = synthetic_rows()
        result = analysis.analyze_endpoint(rows, OUTCOME)
        for row in result["contrasts"]:
            row.update(p_raw=.0001, p_adjusted=.2, stars="****")
        with contextlib.redirect_stdout(io.StringIO()):
            fig = plots.region_figure(rows, [OUTCOME], "%", statistics={OUTCOME: result})
        self.addCleanup(plots.pyplot().close, fig)
        self.assertFalse(any(text.get_text().strip("*") == "" for text in fig.axes[0].texts))
        result["metadata"]["status"] = "insufficient_data"
        fig2 = plots.lesion_figure(rows, "mm3", statistics={"Total stroke volume [mm3]": result})
        self.addCleanup(plots.pyplot().close, fig2)
        self.assertTrue(any("Statistics unavailable" in text.get_text() for text in fig2.texts))
        self.assertFalse(any("No significant comparisons" in text.get_text() for text in fig2.texts))

    def test_significant_brackets_are_above_data_on_both_axis_scales(self):
        rows = synthetic_rows()
        result = analysis.analyze_endpoint(rows, OUTCOME)
        with contextlib.redirect_stdout(io.StringIO()):
            plain = plots.region_figure(rows, [OUTCOME], "%")
            annotated = plots.region_figure(rows, [OUTCOME], "%", statistics={OUTCOME: result})
        self.addCleanup(plots.pyplot().close, plain)
        self.addCleanup(plots.pyplot().close, annotated)
        stars = [text for text in annotated.axes[0].texts if text.get_text().startswith("*")]
        self.assertEqual(len(stars), sum(row["p_adjusted"] < .05 for row in result["contrasts"]))
        self.assertTrue(all(max(row[OUTCOME] for row in rows) < text.get_position()[1] < annotated.axes[0].get_ylim()[1]
                            for text in stars))
        total = analysis.analyze_endpoint(rows, "Total stroke volume [mm3]")
        fig = plots.lesion_figure(rows, "mm3", statistics={"Total stroke volume [mm3]": total})
        self.addCleanup(plots.pyplot().close, fig)
        self.assertTrue(any(text.get_text().startswith("*") for text in fig.axes[0].texts))

    def test_plot_and_statistics_export_include_results_and_allow_missing_cells(self):
        rows = synthetic_rows()
        rows[0][OUTCOME] = None
        path = self.folder / "cohort.csv"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        original = path.read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):
            paths = plots.plot_csv(path, read_input=mock.Mock(side_effect=[OUTCOME, "mm3"]), statistics=True)
        self.assertEqual(len(paths), 7)
        folder = self.folder / "cohort_plots"
        with (folder / "mixed_model_status.csv").open(newline="") as stream:
            status = list(csv.DictReader(stream))
        self.assertEqual(len(status), 2)
        self.assertEqual(status[0]["n_missing_values"], "1")
        self.assertTrue(all(row["status"] in ("ok", "simplified_ols") for row in status), status)
        with (folder / "mixed_model_contrasts.csv").open(newline="") as stream:
            comparisons = list(csv.DictReader(stream))
        self.assertEqual(len(comparisons), 8)
        self.assertTrue(all(row["family_size"] == "4" for row in comparisons))
        self.assertEqual(path.read_bytes(), original)

    def test_statistics_flag_requires_plot_and_group_before_export(self):
        base = ["-i", str(self.folder), "-a", "parental.csv", "-o", str(self.folder / "cohort.csv"), "--statistics"]
        for extra in ([], ["--plot"], ["--group", "groups.csv"]):
            with self.subTest(extra=extra), contextlib.redirect_stderr(io.StringIO()) as stderr:
                with self.assertRaises(SystemExit) as caught:
                    aggregate.main([*base, *extra])
                self.assertEqual(caught.exception.code, 2)
                self.assertIn("requires both --plot and --group", stderr.getvalue())
                self.assertFalse((self.folder / "cohort.csv").exists())


if __name__ == "__main__":
    unittest.main()
