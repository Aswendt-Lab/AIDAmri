"""Exercise acronym report rendering and batch integration with small NIfTI images."""

import contextlib
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))

import batchProc
from helper_tools import batch_qc_reports as reports
from matplotlib.axes import Axes
from matplotlib.figure import Figure


class EntryParser(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.entries = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "div" and attrs.get("class") == "report-entry":
            self.entries.append(attrs)


class AcronymReportTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.project = Path(self.temp_dir.name)
        folder = self.project / "sub-01" / "ses-01" / "anat"
        folder.mkdir(parents=True)
        self.bet = folder / "scanBet.nii.gz"
        self.anno = folder / "scanBet_AnnoSplit_parental.nii.gz"
        nib.save(nib.Nifti1Image(np.ones((4, 4, 4)), np.eye(4)), self.bet)
        # Include both PTLp labels and unrelated atlas labels.
        self.labels = np.broadcast_to(
            np.array([22, 2022, 891, 892], dtype=np.int16), (4, 4, 4)
        ).copy()
        nib.save(nib.Nifti1Image(self.labels, np.eye(4)), self.anno)

    def test_acronym_filters_overlay_and_registration_keeps_all_labels(self):
        imshow = Axes.imshow
        displayed_labels = set()

        def capture_overlay(ax, data, *args, **kwargs):
            if np.ma.isMaskedArray(data):
                displayed_labels.update(np.unique(data.compressed()))
            return imshow(ax, data, *args, **kwargs)

        with patch.object(Axes, "imshow", capture_overlay), \
                patch.object(Figure, "suptitle", autospec=True) as title:
            html_path, count = reports.build_acronym_qc_report(
                self.project, "ptlp", n_slices=1,
                custom_parameters=[("--t2-frac", 0.1)],
            )
        self.assertEqual(title.call_args.args[1], f"PTLp Report: {self.bet.name} + {self.anno.name}")
        self.assertEqual(displayed_labels, {22, 2022})
        self.assertEqual(count, 1)
        self.assertEqual(html_path, self.project / "Report" / "PTLp" / "PTLp_report.html")
        html = html_path.read_text()
        self.assertNotIn("<select id='groupDropdown'", html)
        expected_title = "PTLp Report BET + AnnoSplit_parental (labels 22, 2022)"
        self.assertIn(f"<title>{expected_title}</title>", html)
        self.assertIn(f"<h1>{expected_title}</h1>", html)
        self.assertIn("--t2-frac", html)
        pngs = list(html_path.parent.glob("*.png"))
        self.assertEqual(len(pngs), 1)
        self.assertIn(pngs[0].name, html)
        self.assertTrue(pngs[0].read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

        displayed_labels.clear()
        with patch.object(Axes, "imshow", capture_overlay):
            registration_path, count = reports.build_registration_qc_report(
                self.project, n_slices=1
            )
        self.assertEqual(displayed_labels, {22, 2022, 891, 892})
        self.assertEqual(count, 1)
        self.assertTrue(registration_path.is_file())
        np.testing.assert_array_equal(nib.load(self.anno).get_fdata(), self.labels)

    def test_missing_bet_and_incompatible_shapes_are_skipped(self):
        self.bet.unlink()
        with self.assertLogs(level="WARNING"):
            self.assertEqual(reports.build_acronym_qc_report(self.project, "ptlp"), (None, 0))
        nib.save(nib.Nifti1Image(np.ones((3, 3, 3)), np.eye(4)), self.bet)
        with self.assertLogs(level="WARNING"):
            self.assertEqual(reports.build_acronym_qc_report(self.project, "ptlp"), (None, 0))
        self.assertFalse((self.project / "Report" / "PTLp" / "PTLp_report.html").exists())

    def test_multiple_regions_share_overlay_report_and_canonical_names(self):
        labels = np.broadcast_to(
            np.array([22, 2022, 31, 2031, 891, 892], dtype=np.int16), (4, 4, 6)
        ).copy()
        nib.save(nib.Nifti1Image(np.ones(labels.shape), np.eye(4)), self.bet)
        nib.save(nib.Nifti1Image(labels, np.eye(4)), self.anno)
        imshow = Axes.imshow
        displayed_labels = set()

        def capture_overlay(ax, data, *args, **kwargs):
            if np.ma.isMaskedArray(data):
                displayed_labels.update(np.unique(data.compressed()))
            return imshow(ax, data, *args, **kwargs)

        with patch.object(Axes, "imshow", capture_overlay), \
                patch.object(Figure, "suptitle", autospec=True) as title:
            html_path, count = reports.build_acronym_qc_report(
                self.project, ["ptlp", "aca", "PTLP"], n_slices=1
            )
        self.assertEqual(displayed_labels, {22, 2022, 31, 2031})
        self.assertEqual(count, 1)
        self.assertEqual(html_path, self.project / "Report" / "PTLp_ACA" / "PTLp_ACA_report.html")
        self.assertEqual(title.call_args.args[1], f"PTLp + ACA Report: {self.bet.name} + {self.anno.name}")
        expected_title = "PTLp + ACA Report BET + AnnoSplit_parental (labels 22, 2022, 31, 2031)"
        html = html_path.read_text()
        self.assertIn(f"<title>{expected_title}</title>", html)
        self.assertIn(f"<h1>{expected_title}</h1>", html)
        self.assertEqual(len(list(html_path.parent.glob("*.png"))), 1)
        np.testing.assert_array_equal(nib.load(self.anno).get_fdata(), labels)

    def test_cli_accepts_multiple_values_and_repeated_acronym_option(self):
        for selection in ("acronym", "all"):
            for options in (
                ["--acronym", "ptlp", "ACA", "PTLP"],
                ["--acronym", "ptlp", "--acronym", "ACA", "PTLP"],
            ):
                with self.subTest(selection=selection, options=options), \
                        contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(reports.main([
                        "-i", str(self.project), "--report", selection, "--n-slices", "1", *options,
                    ]), 0)
                    created = {p.parent.name for p in (self.project / "Report").glob("*/*.html")}
                    expected = {"PTLp_ACA"}
                    if selection == "all":
                        expected.update({"BET", "Registration"})
                    self.assertEqual(created, expected)

    def test_cli_acronym_and_all_generate_expected_reports(self):
        for selection, expected in (("acronym", {"PTLp"}), ("all", {"BET", "Registration", "PTLp"})):
            with self.subTest(selection=selection), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(reports.main([
                    "-i", str(self.project), "--report", selection, "--n-slices", "1",
                    "--acronym", "pTlP",
                ]), 0)
                created = {p.parent.name for p in (self.project / "Report").glob("*/*.html")}
                self.assertEqual(created, expected)

    def test_lookup_preserves_file_spelling_and_adds_2000(self):
        for query, expected in (
            ("ptlp", ("PTLp", (22, 2022))),
            ("PTLP", ("PTLp", (22, 2022))),
            ("CC", ("cc", (776, 2776))),
            ("ssp-bfd", ("SSp-bfd", (329, 2329))),
        ):
            with self.subTest(query=query):
                self.assertEqual(reports._resolve_acronym(query), expected)

    def test_cli_defaults_to_cc_for_all_and_acronym_reports(self):
        for options, expected in (
            (["--report", "acronym"], {"cc"}),
            ([], {"BET", "Registration", "cc"}),
            (["--report", "all"], {"BET", "Registration", "cc"}),
        ):
            with self.subTest(options=options), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(reports.main([
                    "-i", str(self.project), "--n-slices", "1", *options,
                ]), 0)
            created = {path.parent.name for path in (self.project / "Report").glob("*/*.html")}
            self.assertEqual(created, expected)
            html_path = self.project / "Report" / "cc" / "cc_report.html"
            self.assertIn("cc Report BET + AnnoSplit_parental (labels 776, 2776)", html_path.read_text())

    def test_cli_rejects_empty_or_unknown_acronym_before_writing_reports(self):
        for options, message in (
            (["--acronym", "unknown-region"], "Unknown acronym"),
            (["--acronym", "ptlp", "unknown-region"], "Unknown acronym"),
            (["--acronym"], "expected at least one argument"),
        ):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()) as stderr:
                with self.assertRaises(SystemExit) as error:
                    reports.main(["-i", str(self.project), *options])
                self.assertEqual(error.exception.code, 2)
                self.assertIn(message, stderr.getvalue())
                self.assertFalse((self.project / "Report").exists())

    def test_multiple_acronym_lookup_preserves_order_and_removes_duplicates(self):
        self.assertEqual(
            reports._resolve_acronyms(["ptlp", "ACA", "PTLP", "hip"]),
            (["PTLp", "ACA", "HIP"], (22, 2022, 31, 2031, 1080, 3080)),
        )
        with self.assertRaisesRegex(ValueError, "At least one"):
            reports._resolve_acronyms([])

    def test_atlas_selects_lookup_table_and_parental_is_default(self):
        self.assertEqual(
            reports._resolve_acronyms(["mop2/3", "FRP5"], atlas="detailed"),
            (["MOp2/3", "FRP5"], (943, 2943, 526157192, 526159192)),
        )
        with self.assertRaisesRegex(ValueError, "acronym_rsfMRI.txt"):
            reports._resolve_acronym("mop2/3")
        self.assertEqual(
            reports._resolve_acronym("ptlp"),
            reports._resolve_acronym("ptlp", atlas="parental"),
        )
        args = reports._build_argument_parser().parse_args(["-i", str(self.project)])
        self.assertEqual(args.atlas, "parental")

    def test_detailed_lookup_handles_acronyms_containing_spaces(self):
        # These rows precede int in the real table; unrestricted split breaks
        # the lookup before int can be reached.
        for query, expected in (
            ("INT", ("int", (6, 2006))),
            ("cul4, 5", ("CUL4, 5", (1091, 3091))),
            ("FIBER TRACTS", ("fiber tracts", (1009, 3009))),
        ):
            with self.subTest(query=query):
                self.assertEqual(reports._resolve_acronym(query, atlas="detailed"), expected)

    def test_exact_lookup_resolves_every_entry_in_both_atlas_tables(self):
        for atlas, (filename, _) in reports.ATLAS_CONFIG.items():
            source = Path(reports.__file__).resolve().parents[2] / "lib" / filename
            for line in source.read_text(encoding="utf-8").splitlines():
                if not line.strip():
                    continue
                label, acronym = line.strip().split(maxsplit=1)
                with self.subTest(atlas=atlas, acronym=acronym):
                    self.assertEqual(
                        reports._resolve_acronym(acronym, atlas=atlas),
                        (acronym, (int(label), int(label) + 2000)),
                    )

    def test_case_distinct_regions_can_be_selected_together(self):
        self.assertEqual(
            reports._resolve_acronyms(["CM", "cm", "ipf", "IPF", "CM"], atlas="detailed"),
            (["CM", "cm", "ipf", "IPF"], (599, 2599, 967, 2967, 49, 2049, 624, 2624)),
        )
        self.assertEqual(reports._resolve_acronym("CM"), ("cm", (967, 2967)))

    def test_cli_rejects_ambiguous_case_with_candidates_before_creating_reports(self):
        for query, candidates in (
            ("cM", ("'CM' (label 599)", "'cm' (label 967)")),
            ("iPf", ("'ipf' (label 49)", "'IPF' (label 624)")),
        ):
            with self.subTest(query=query), contextlib.redirect_stderr(io.StringIO()) as stderr:
                with self.assertRaises(SystemExit) as error:
                    reports.main([
                        "-i", str(self.project), "--report", "all",
                        "--atlas", "detailed", "--acronym", "ptlp", query,
                    ])
                self.assertEqual(error.exception.code, 2)
                self.assertIn("Ambiguous acronym", stderr.getvalue())
                for candidate in candidates:
                    self.assertIn(candidate, stderr.getvalue())
                self.assertIn("Use the exact spelling", stderr.getvalue())
                self.assertFalse((self.project / "Report").exists())

    def test_detailed_cli_uses_detailed_file_and_renders_multiple_regions(self):
        detailed_path = self.anno.with_name("scanBet_AnnoSplit.nii.gz")
        labels = np.broadcast_to(
            np.array([943, 2943, 68, 2068, 22, 2022], dtype=np.int16), (4, 4, 6)
        ).copy()
        nib.save(nib.Nifti1Image(np.ones(labels.shape), np.eye(4)), self.bet)
        nib.save(nib.Nifti1Image(labels, np.eye(4)), detailed_path)
        imshow = Axes.imshow
        displayed_labels = set()

        def capture_overlay(ax, data, *args, **kwargs):
            if np.ma.isMaskedArray(data):
                displayed_labels.update(np.unique(data.compressed()))
            return imshow(ax, data, *args, **kwargs)

        with patch.object(Axes, "imshow", capture_overlay), \
                patch.object(Figure, "suptitle", autospec=True) as title, \
                contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(reports.main([
                "-i", str(self.project), "--report", "acronym", "--atlas", "detailed",
                "--acronym", "mop2/3", "frp1", "--n-slices", "1",
            ]), 0)
        self.assertEqual(displayed_labels, {943, 2943, 68, 2068})
        self.assertEqual(title.call_args.args[1], f"MOp2/3 + FRP1 Report: {self.bet.name} + {detailed_path.name}")
        html_path = self.project / "Report" / "MOp2_3_FRP1" / "MOp2_3_FRP1_detailed_report.html"
        html = html_path.read_text()
        self.assertIn("MOp2/3 + FRP1 Report BET + AnnoSplit (labels 943, 2943, 68, 2068)", html)
        self.assertIn(detailed_path.name, html)
        self.assertNotIn(self.anno.name, html)
        pngs = list(html_path.parent.glob("*.png"))
        self.assertEqual(len(pngs), 1)
        self.assertIn(pngs[0].name, html)
        np.testing.assert_array_equal(nib.load(detailed_path).get_fdata(), labels)

    def test_detailed_does_not_fall_back_to_parental_annotation(self):
        self.assertEqual(
            reports.build_acronym_qc_report(self.project, "ptlp", atlas="detailed"),
            (None, 0),
        )

    def test_detailed_all_preserves_parental_reports(self):
        parental_html, _ = reports.build_acronym_qc_report(
            self.project, "ptlp", n_slices=1, atlas="parental"
        )
        parental_files = {p: p.read_bytes() for p in parental_html.parent.iterdir()}
        detailed_path = self.anno.with_name("scanBet_AnnoSplit.nii.gz")
        nib.save(nib.Nifti1Image(self.labels, np.eye(4)), detailed_path)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(reports.main([
                "-i", str(self.project), "--report", "all", "--atlas", "detailed",
                "--acronym", "ptlp", "--n-slices", "1",
            ]), 0)
        for path, content in parental_files.items():
            self.assertEqual(path.read_bytes(), content)
        detailed_html = parental_html.with_name("PTLp_detailed_report.html")
        self.assertIn("PTLp Report BET + AnnoSplit (labels 22, 2022)", detailed_html.read_text())
        registration_html = self.project / "Report" / "Registration" / "registration_report.html"
        self.assertIn(self.anno.name, registration_html.read_text())
        self.assertTrue((self.project / "Report" / "BET" / "bet_report.html").is_file())

    def test_bet_and_registration_cli_do_not_require_acronym(self):
        for selection, builder in (
            ("bet", "build_bet_qc_report"),
            ("registration", "build_registration_qc_report"),
        ):
            with self.subTest(selection=selection), \
                    patch.object(reports, builder, return_value=(None, 0)) as build, \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(reports.main(["-i", str(self.project), "--report", selection]), 0)
                build.assert_called_once()

    def test_batch_compatibility_wrapper_uses_acronym_file(self):
        with patch.object(reports, "build_acronym_qc_report", return_value=(None, 0)) as build:
            self.assertEqual(reports.build_cc_qc_report(self.project, n_slices=7), (None, 0))
            build.assert_called_once_with(self.project, "cc", n_slices=7, custom_parameters=None, group=None)

    def test_cli_groups_apply_to_all_qc_reports_with_unassigned_subjects_kept(self):
        for subject in ("sub-02", "sub-03"):
            folder = self.project / subject / "ses-01" / "anat"
            folder.mkdir(parents=True)
            nib.save(nib.Nifti1Image(np.ones((4, 4, 4)), np.eye(4)), folder / self.bet.name)
            nib.save(nib.Nifti1Image(self.labels, np.eye(4)), folder / self.anno.name)
        csv_path = self.project / "groups.csv"
        csv_path.write_text("tDCS,Sham\nsub-01,sub-02\n", encoding="utf-8")
        for flag in ("-g", "--group"):
            with self.subTest(flag=flag), contextlib.redirect_stdout(io.StringIO()):
                status = reports.main([
                    "-i", str(self.project), "--report", "all", "--acronym", "ptlp",
                    "--n-slices", "1", flag, str(csv_path),
                ])
            self.assertEqual(status, 0)
            html_files = list((self.project / "Report").glob("*/*.html"))
            self.assertEqual({path.parent.name for path in html_files}, {"BET", "Registration", "PTLp"})
            for path in html_files:
                html = path.read_text()
                self.assertIn("<select id='groupDropdown'", html)
                self.assertIn("<option value='group-0'>tDCS</option>", html)
                self.assertIn("<option value='group-1'>Sham</option>", html)
                self.assertIn("<option value='unassigned'>Unassigned</option>", html)
                memberships = {entry["data-subject"]: json.loads(entry["data-groups"])
                               for entry in EntryParser(html).entries}
                self.assertEqual(memberships, {"sub-01": ["group-0"], "sub-02": ["group-1"],
                                               "sub-03": ["unassigned"]})

    def test_python_report_builders_accept_group_csv_including_cc_and_detailed_atlas(self):
        csv_path = self.project / "groups.csv"
        csv_path.write_text("tDCS,Sham\nsub-01,\n", encoding="utf-8")
        detailed_anno = self.anno.with_name("scanBet_AnnoSplit.nii.gz")
        nib.save(nib.Nifti1Image(self.labels, np.eye(4)), detailed_anno)
        for builder, kwargs in (
            (reports.build_bet_qc_report, {}),
            (reports.build_registration_qc_report, {}),
            (reports.build_acronym_qc_report, {"acronym": "ptlp", "atlas": "detailed"}),
            (reports.build_cc_qc_report, {}),
        ):
            with self.subTest(builder=builder.__name__):
                html_path, count = builder(self.project, n_slices=1, group=csv_path, **kwargs)
            self.assertEqual(count, 1)
            html = html_path.read_text()
            self.assertIn("<select id='groupDropdown'", html)
            self.assertEqual(json.loads(EntryParser(html).entries[0]["data-groups"]), ["group-0"])

    def test_invalid_group_csv_is_rejected_before_any_report_is_created(self):
        csv_path = self.project / "groups.csv"
        for content in (None, "tDCS,tDCS\n", 'tDCS,Sham\n"unterminated'):
            if content is not None:
                csv_path.write_text(content, encoding="utf-8")
            with self.subTest(content=content), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    reports.main([
                        "-i", str(self.project), "--report", "all", "--acronym", "ptlp",
                        "-g", str(csv_path),
                    ])
                self.assertEqual(error.exception.code, 2)
            self.assertFalse((self.project / "Report").exists())

    def test_batch_creates_cc_only_for_registration_and_survives_report_failure(self):
        parameters = [("--t2-frac", 0.1)]
        with patch.object(reports, "build_cc_qc_report", return_value=(None, 0)) as cc, \
                patch.object(reports, "build_bet_qc_report", return_value=(None, 0)), \
                patch.object(reports, "build_registration_qc_report", side_effect=RuntimeError("failed")), \
                contextlib.redirect_stdout(io.StringIO()):
            batchProc.create_qc_reports(self.project, ["preprocess", "process"], parameters)
            cc.assert_not_called()
            with self.assertLogs(level="WARNING"):
                batchProc.create_qc_reports(self.project, ["registration"], parameters)
            cc.assert_called_once_with(self.project, n_slices=7, custom_parameters=parameters)


if __name__ == "__main__":
    unittest.main()
