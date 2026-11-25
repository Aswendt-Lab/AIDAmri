"""Check interactive selection, plotted statistics and headless PNG/PDF export."""

import contextlib
import csv
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]


def load_script(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin/3.1_T2Processing" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


plots = load_script("affected_regions_plots", "affected_regions_plots.py")
aggregate = load_script("plotting_cohort_export", "aggregate_affected_regions.py")


class AffectedRegionPlotTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)
        self.csv = self.folder / "cohort.csv"
        self.columns = ["L_A [mm3]", "L_B [mm3]", "L_A [%]", "L_B [%]"]
        self.rows = []
        for index, (session, group, value) in enumerate([
                ("ses-PT3", "tDCS", 0), ("ses-PT3", "tDCS", 1),
                ("ses-PT3", "tDCS", 3), ("ses-PT14", "Sham", 5)]):
            self.rows.append({"subject": f"sub-{index}", "session": session, "group": group,
                              plots.TOTAL_COLUMNS["mm3"]: value, plots.TOTAL_COLUMNS["%"]: value / 2,
                              "L_A [mm3]": value, "L_B [mm3]": 0,
                              "L_A [%]": value / 2, "L_B [%]": 0})

    def write_csv(self, rows=None, with_group=True):
        header = ["subject", "session"] + (["group"] if with_group else [])
        header += [*plots.TOTAL_COLUMNS.values(), *self.columns]
        with self.csv.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=header, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(self.rows if rows is None else rows)

    def source_folder(self):
        root = self.folder / "proc"
        path = root / "sub-A/ses-PT3/anat/affected_Regions/scan_affectedRegions_Parental.csv"
        path.parent.mkdir(parents=True)
        path.write_text("Label_id,Brain region,Affected region percentage [%],Region stroke volume [mm3],,"
                        "Stroke percent of brain [%],Total stroke volume [mm3]\n"
                        "22,L_Region,2.50,0.1000,,0.50,1.0000\n", encoding="utf-8")
        return root

    def test_explicit_selection_preserves_order_deduplicates_and_handles_commas_in_names(self):
        long_name = "L_Primary_somatosensory_area,_lower_limb [%]"
        columns = [*self.columns, long_name]
        reader = mock.Mock(return_value=f"{long_name}, L_B [mm3], {long_name}, L_A [%]")
        with contextlib.redirect_stdout(io.StringIO()):
            result = plots.select_regions(columns, reader)
        self.assertEqual(result, [long_name, "L_B [mm3]", "L_A [%]"])
        reader.assert_called_once()

    def test_all_requests_unit_and_uses_only_regional_columns_of_that_unit(self):
        reader = mock.Mock(side_effect=["all", "invalid", "%"])
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            result = plots.select_regions(self.columns, reader)
        self.assertEqual(result, ["L_A [%]", "L_B [%]"])
        self.assertEqual(reader.call_count, 3)
        self.assertIn("Please enter mm3 or %", stdout.getvalue())

    def test_invalid_explicit_selections_are_requested_again(self):
        reader = mock.Mock(side_effect=["", "L_missing [%]", "L_A [%],", "L_A [%]"])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(plots.select_regions(self.columns, reader), ["L_A [%]"])
        self.assertEqual(reader.call_count, 4)

    def test_read_csv_normalizes_missing_groups_and_keeps_optional_group_absent(self):
        self.rows[0]["group"] = ""
        self.write_csv()
        columns, rows = plots.read_cohort_csv(self.csv)
        self.assertEqual(columns, self.columns)
        self.assertEqual(rows[0]["group"], "Unassigned")
        self.assertEqual(plots.group_layout(rows)[0], ["ses-PT3", "ses-PT14"])
        self.write_csv(with_group=False)
        _, rows = plots.read_cohort_csv(self.csv)
        self.assertEqual(plots.group_layout(rows), (["ses-PT3", "ses-PT14"], [""]))

    def test_invalid_numeric_csv_data_is_rejected_before_plotting(self):
        for value in ("NaN", "inf", "", "bad", -1):
            with self.subTest(value=value):
                self.rows[0]["L_A [%]"] = value
                self.write_csv()
                with self.assertRaises(ValueError) as caught:
                    plots.read_cohort_csv(self.csv)
                self.assertIn(str(self.csv), str(caught.exception))

    def test_box_statistics_use_min_max_and_only_regional_statistics_drop_zero(self):
        self.assertIsNone(plots.box_statistics([0, 0], positive_only=True))
        regional = plots.box_statistics([0, 1, 2, 3, 4, 100], positive_only=True)
        self.assertEqual((regional["q1"], regional["med"], regional["q3"]), (2, 3, 4))
        self.assertEqual((regional["whislo"], regional["whishi"]), (1, 100))
        lesion = plots.box_statistics([0, 1, 2, 3, 4])
        self.assertEqual((lesion["q1"], lesion["med"], lesion["q3"]), (1, 2, 3))
        self.assertEqual((lesion["whislo"], lesion["whishi"]), (0, 4))

    def test_regional_plot_groups_sessions_labels_units_and_excluded_zeros(self):
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout), mock.patch.object(plots, "draw_box", wraps=plots.draw_box) as draw:
            fig = plots.region_figure(self.rows, ["L_A [%]", "L_B [%]"], "%")
        self.addCleanup(plots.pyplot().close, fig)
        ax = fig.axes[0]
        self.assertEqual(ax.get_yscale(), "log")
        self.assertEqual(ax.get_ylabel(), "% infarct area")
        self.assertEqual([label.get_text() for label in ax.get_xticklabels()], ["L_A", "L_B"])
        self.assertEqual([text.get_text() for text in ax.get_legend().get_texts()],
                         ["tDCS PT3", "Sham PT3", "tDCS PT14", "Sham PT14"])
        stats = [call.args[1] for call in draw.call_args_list]
        self.assertEqual([(item["whislo"], item["whishi"]) for item in stats], [(.5, 1.5), (2.5, 2.5)])
        self.assertIn("excluded zeros: 5", fig.texts[0].get_text())
        self.assertTrue(any(text.get_text() == "No positive\nvalues" for text in ax.texts))
        self.assertEqual(len(ax.collections), 1)  # The single-observation box keeps a colored point.
        self.assertIn("L_A [%] / tDCS / ses-PT3: 2 positive values; excluded zeros: 1", stdout.getvalue())
        self.assertNotEqual(ax.get_legend().legend_handles[0].get_facecolor(),
                            ax.get_legend().legend_handles[2].get_facecolor())

    def test_entirely_zero_regional_plot_is_skipped(self):
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            fig = plots.region_figure(self.rows, ["L_B [%]"], "%")
        self.assertIsNone(fig)
        self.assertIn("Skipping % regional plot: no positive values (excluded zeros: 4)", stdout.getvalue())

    def test_lesion_plot_includes_zero_and_all_points_with_reproducible_jitter(self):
        with mock.patch.object(plots, "draw_box", wraps=plots.draw_box) as draw:
            fig = plots.lesion_figure(self.rows, "mm3")
        self.addCleanup(plots.pyplot().close, fig)
        ax = fig.axes[0]
        self.assertEqual(ax.get_yscale(), "linear")
        self.assertEqual(ax.get_xlabel(), "Time point")
        self.assertEqual(ax.get_ylabel(), "Stroke volume (mm³)")
        self.assertEqual([label.get_text() for label in ax.get_xticklabels()], ["PT3", "PT14"])
        stats = [call.args[1] for call in draw.call_args_list]
        self.assertEqual((stats[0]["whislo"], stats[0]["med"], stats[0]["whishi"]), (0, 1, 3))
        points = [collection.get_offsets().tolist() for collection in ax.collections]
        self.assertEqual(sum(len(values) for values in points), 4)
        second = plots.lesion_figure(self.rows, "mm3")
        self.addCleanup(plots.pyplot().close, second)
        self.assertEqual(points, [collection.get_offsets().tolist() for collection in second.axes[0].collections])

    def test_plot_export_writes_png_pdf_and_prints_statistics_definition(self):
        self.write_csv()
        original = self.csv.read_bytes()
        stdout = io.StringIO()
        reader = mock.Mock(side_effect=["all", "%", "mm3"])
        with contextlib.redirect_stdout(stdout):
            paths = plots.plot_csv(self.csv, reader)
        self.assertEqual(reader.call_count, 3)
        self.assertEqual({path.name for path in paths},
                         {"regions_percent.png", "regions_percent.pdf", "lesion_size_mm3.png", "lesion_size_mm3.pdf"})
        for path in paths:
            self.assertEqual(path.parent, self.folder / "cohort_plots")
            self.assertTrue(path.read_bytes().startswith(b"%PDF" if path.suffix == ".pdf" else b"\x89PNG"))
        self.assertIn("whiskers = minimum to maximum", stdout.getvalue())
        self.assertEqual(self.csv.read_bytes(), original)

    def test_mixed_units_generate_two_regions_plots_without_additional_unit_prompt(self):
        self.write_csv(with_group=False)
        reader = mock.Mock(side_effect=["L_A [%],L_A [mm3]", "%"])
        with contextlib.redirect_stdout(io.StringIO()), mock.patch.object(plots, "save_figure", return_value=[]) as save:
            plots.plot_csv(self.csv, reader)
        self.assertEqual(reader.call_count, 2)
        self.assertEqual([call.args[2] for call in save.call_args_list],
                         ["regions_mm3", "regions_percent", "lesion_size_percent"])

    def test_wide_png_resolution_is_limited_and_vector_pdf_is_kept(self):
        fig = mock.Mock()
        fig.get_size_inches.return_value = (1000, 6)
        fig.savefig.side_effect = lambda path, **kwargs: path.write_bytes(b"figure")
        with contextlib.redirect_stdout(io.StringIO()):
            paths = plots.save_figure(fig, self.folder, "wide")
        self.assertEqual([path.suffix for path in paths], [".png", ".pdf"])
        self.assertLess(fig.savefig.call_args_list[0].kwargs["dpi"], 300)

    def test_cli_plot_runs_prompts_and_saves_images(self):
        root = self.source_folder()
        with mock.patch("builtins.input", side_effect=["all", "%", "mm3"]), contextlib.redirect_stdout(io.StringIO()):
            aggregate.main(["-i", str(root), "-a", "parental.csv", "-o", str(self.csv), "--plot"])
        self.assertTrue(self.csv.exists())
        self.assertTrue((self.folder / "cohort_plots/regions_percent.png").exists())
        self.assertTrue((self.folder / "cohort_plots/lesion_size_mm3.pdf").exists())

    def test_no_plot_does_not_prompt_or_import_plotting(self):
        root = self.source_folder()
        with mock.patch("builtins.input", side_effect=AssertionError("Unexpected input")), \
                mock.patch("builtins.__import__", wraps=__import__) as imports, \
                contextlib.redirect_stdout(io.StringIO()):
            aggregate.main(["-i", str(root), "-a", "parental.csv", "-o", str(self.csv)])
        self.assertFalse(any(call.args[0].startswith(("matplotlib", "affected_regions_plots"))
                             for call in imports.call_args_list))

    def test_plot_failure_or_cancel_keeps_exported_csv(self):
        root = self.source_folder()
        for error in (EOFError(), KeyboardInterrupt(), RuntimeError("plot failed"), ImportError("no Matplotlib")):
            with self.subTest(error=type(error).__name__), \
                    mock.patch.object(plots, "plot_csv", side_effect=error), \
                    contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as stderr:
                with self.assertRaises(SystemExit) as caught:
                    aggregate.main(["-i", str(root), "-a", "parental.csv", "-o", str(self.csv), "--plot"])
                self.assertEqual(caught.exception.code, 2)
                self.assertIn("CSV retained:", stderr.getvalue())
                self.assertIn("0.1000", self.csv.read_text())


if __name__ == "__main__":
    unittest.main()
