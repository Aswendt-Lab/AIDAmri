"""Check suffix discovery, overlay safety and rendering using small NIfTI images."""

import contextlib
import io
import json
from html.parser import HTMLParser
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from helper_tools import display_nifti_report as reports
from matplotlib.axes import Axes
from matplotlib.figure import Figure


class ReportHtml(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.options = {}
        self.entries = []
        self.select = None
        self.option = None
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "select":
            self.select = attrs["id"]
            self.options[self.select] = []
        elif tag == "option":
            self.option = [attrs["value"], ""]
        elif tag == "div" and attrs.get("class") == "report-entry":
            self.entries.append(attrs)

    def handle_data(self, data):
        if self.option is not None:
            self.option[1] += data

    def handle_endtag(self, tag):
        if tag == "option":
            self.options[self.select].append(tuple(self.option))
            self.option = None
        elif tag == "select":
            self.select = None


class DisplayNiftiReportTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.project = Path(self.temp_dir.name)
        self.modality = self.project / "sub-NR1037" / "ses-PT3" / "dwi"
        self.base = self.save(self.modality / "DSI_studio" / "scan.fz.fa.nii.gz")

    def save(self, path, data=None, affine=None):
        path.parent.mkdir(parents=True, exist_ok=True)
        if data is None:
            data = np.arange(4 * 5 * 6, dtype=np.float32).reshape(4, 5, 6)
        image = nib.Nifti1Image(data, np.eye(4) if affine is None else affine)
        nib.save(image, path)
        return path

    def overlay(self, folder=None, affine=None):
        labels = np.zeros((4, 5, 6), dtype=np.int16)
        labels[1:3, 1:4, 1:5] = 22
        labels[2, 2, 2] = 2022
        return self.save(
            (folder or self.modality) / "scan_AnnoSplit_parental.nii.gz",
            data=labels, affine=affine,
        )

    def build(self, **kwargs):
        return reports.build_display_nifti_report(
            self.project, "fz.fa.nii.gz", n_slices=1, **kwargs,
        )

    def test_single_input_searches_all_subjects_sessions_modalities_and_subfolders(self):
        self.save(self.modality / "other" / "second.fz.fa.nii.gz")
        self.save(self.project / "sub-02" / "ses-02" / "custom" / "scan.fz.fa.nii.gz")
        self.save(self.project / "sub-NR1037" / "ses-PT4" / "anat" / "scan.fz.fa.nii.gz")
        # Files outside Subject/Session/Modality and other suffixes are excluded.
        self.save(self.project / "Report" / "scan.fz.fa.nii.gz")
        self.save(self.modality / "scan.fz.md.nii.gz")
        with patch.object(Figure, "suptitle", autospec=True, side_effect=Figure.suptitle) as title:
            html_path, count = self.build()
        self.assertEqual(len(title.call_args_list), 4)
        self.assertTrue(all(call.args[1].startswith("Display Report: ")
                            for call in title.call_args_list))
        self.assertEqual(count, 4)
        self.assertEqual(html_path, self.project / "Report" / "Display" /
                         "display_nifti_report_fz.fa.nii.gz.html")
        html = html_path.read_text()
        self.assertNotIn("groupDropdown", ReportHtml(html).options)
        self.assertIn("<title>Display Report: fz.fa.nii.gz</title>", html)
        self.assertIn("<h1>Display Report: fz.fa.nii.gz</h1>", html)
        for value in ("subjectDropdown", "sessionDropdown", "modalityDropdown",
                      'id="image-viewer"', 'id="viewer-zoom-in"', "changeZoom",
                      "pointerdown", "wheel", "Created:",
                      "data-modality='custom'", "DSI_studio/scan.fz.fa.nii.gz"):
            self.assertIn(value, html)
        pngs = list(html_path.parent.glob("*.png"))
        self.assertEqual(len(pngs), count)
        for png in pngs:
            self.assertIn(png.name, html)
            self.assertTrue(png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

    def test_modality_selection_only_searches_selected_folders(self):
        for modality in ("anat", "func", "t2map"):
            self.save(self.modality.parent / modality / "nested" / "scan.fz.fa.nii.gz")
        self.save(self.project / "sub-02" / "ses-02" / "dwi" / "nested" / "scan.fz.fa.nii.gz")
        self.save(self.project / "sub-02" / "ses-03" / "dwi" / "scan.fz.fa.nii.gz")
        original_rglob = Path.rglob
        for selection, expected_count in (("all", 6), ("anat", 1), ("dwi", 3),
                                          ("func", 1), ("t2map", 1)):
            searched = []

            def capture_rglob(path, pattern):
                searched.append(path.name)
                return original_rglob(path, pattern)

            with self.subTest(selection=selection), patch.object(Path, "rglob", capture_rglob):
                html_path, count = self.build(modality=selection)
            self.assertEqual(count, expected_count)
            if selection != "all":
                self.assertEqual(set(searched), {selection})
                html = html_path.read_text()
                for candidate in ("anat", "dwi", "func", "t2map"):
                    self.assertEqual(f"data-modality='{candidate}'" in html, candidate == selection)

    def test_overlay_pairs_subfolders_and_uses_requested_opacity(self):
        overlay = self.overlay()
        images = []
        original_imshow = Axes.imshow

        def capture(ax, data, *args, **kwargs):
            images.append((data.copy(), kwargs))
            return original_imshow(ax, data, *args, **kwargs)

        original = nib.load(overlay).get_fdata()
        with patch.object(Axes, "imshow", capture), \
                patch.object(Figure, "suptitle", autospec=True, side_effect=Figure.suptitle) as title:
            html_path, count = self.build(nifti_in2="AnnoSplit_parental.nii.gz", opacity=60)
        self.assertEqual(title.call_args.args[1], f"Overlay Report: {self.base.name} + {overlay.name}")
        self.assertEqual(count, 1)
        self.assertEqual(html_path, self.project / "Report" / "Overlay" /
                         "overlay_nifti_report_fz.fa.nii.gz_AnnoSplit_parental.nii.gz.html")
        pngs = list(html_path.parent.glob("*.png"))
        self.assertEqual(len(pngs), 1)
        self.assertIn(pngs[0].name, html_path.read_text())
        self.assertEqual(len(images), 6)
        for data, options in images[1::2]:
            self.assertTrue(np.ma.isMaskedArray(data))
            self.assertEqual(options["alpha"], 0.6)
            self.assertEqual(options["cmap"], "tab20")
            self.assertNotIn(0, data.compressed())
        html = html_path.read_text()
        expected_title = "Overlay Report: fz.fa.nii.gz + AnnoSplit_parental.nii.gz"
        self.assertIn(f"<title>{expected_title}</title>", html)
        self.assertIn(f"<h1>{expected_title}</h1>", html)
        self.assertIn("Overlay opacity:</b> 60%", html)
        self.assertIn(overlay.name, html)
        np.testing.assert_array_equal(nib.load(overlay).get_fdata(), original)

    def test_opacity_zero_hides_overlay_and_contours_and_hundred_is_opaque(self):
        self.overlay()
        original_imshow = Axes.imshow
        for opacity in (0, 100):
            with self.subTest(opacity=opacity), \
                    patch.object(Axes, "imshow", autospec=True, side_effect=original_imshow) as imshow, \
                    patch.object(Axes, "contour", autospec=True) as contour:
                _, count = self.build(nifti_in2="AnnoSplit_parental.nii.gz", opacity=opacity)
                self.assertEqual(count, 1)
                overlays = [call for call in imshow.call_args_list if "alpha" in call.kwargs]
                if opacity == 0:
                    self.assertFalse(overlays)
                    contour.assert_not_called()
                else:
                    self.assertEqual(len(overlays), 3)
                    self.assertTrue(all(call.kwargs["alpha"] == 1 for call in overlays))

    def test_different_suffix_selections_preserve_previous_reports_and_images(self):
        self.overlay()
        self.save(self.modality / "scan_OtherAtlas.nii.gz", data=np.ones((4, 5, 6)))
        previous_files = {}
        for nifti_in1, nifti_in2 in (
            ("fz.fa.nii.gz", "AnnoSplit_parental.nii.gz"),
            ("fz.fa.nii.gz", "OtherAtlas.nii.gz"),
            ("fa.nii.gz", "OtherAtlas.nii.gz"),
        ):
            html_path, count = reports.build_display_nifti_report(
                self.project, nifti_in1, nifti_in2=nifti_in2, n_slices=1,
            )
            self.assertEqual(count, 1)
            for path, content in previous_files.items():
                self.assertEqual(path.read_bytes(), content)
            for path in html_path.parent.iterdir():
                previous_files[path] = path.read_bytes()
            html = html_path.read_text()
            matching_pngs = [path for path in previous_files
                             if path.suffix == ".png" and path.name in html]
            self.assertEqual(len(matching_pngs), 1)
        self.assertEqual(len(list(html_path.parent.glob("*.html"))), 3)
        self.assertEqual(len(list(html_path.parent.glob("*.png"))), 3)

    def test_incomplete_pairs_are_not_combined_across_subject_session_or_modality(self):
        for folder in (
            self.project / "sub-02" / "ses-PT3" / "dwi",
            self.project / "sub-NR1037" / "ses-PT4" / "dwi",
            self.project / "sub-NR1037" / "ses-PT3" / "anat",
        ):
            self.overlay(folder=folder)
        with self.assertLogs(level="WARNING") as logs:
            self.assertEqual(self.build(nifti_in2="AnnoSplit_parental.nii.gz"), (None, 0))
        self.assertEqual(len(logs.output), 4)
        self.assertTrue(all("incomplete combination" in line for line in logs.output))
        self.assertFalse(list((self.project / "Report" / "Overlay").glob("*.png")))

    def test_ambiguous_matches_in_either_input_are_warned_and_skipped(self):
        self.overlay()
        for filename in ("duplicate.fz.fa.nii.gz", "duplicate_AnnoSplit_parental.nii.gz"):
            duplicate = self.save(self.modality / "nested" / filename)
            with self.subTest(filename=filename), self.assertLogs(level="WARNING") as logs:
                self.assertEqual(self.build(nifti_in2="AnnoSplit_parental.nii.gz"), (None, 0))
            self.assertIn("ambiguous combination", logs.output[0])
            self.assertIn(filename, logs.output[0])
            duplicate.unlink()

    def test_geometry_mismatches_are_skipped_with_explanation(self):
        for change, expected in (("shape", "Shape mismatch"), ("spacing", "Voxel size mismatch"),
                                 ("translation", "Spatial affine mismatch"),
                                 ("orientation", "Orientation mismatch"),
                                 ("units", "Spatial unit mismatch")):
            affine = np.eye(4)
            if change == "spacing":
                affine[0, 0] = 1.2
            elif change == "translation":
                affine[0, 3] = 0.2
            elif change == "orientation":
                affine[0, 0] = -1
            overlay = self.overlay(affine=affine)
            if change == "shape":
                self.save(overlay, data=np.ones((3, 5, 6), dtype=np.int16))
            elif change == "units":
                base_image = nib.load(self.base)
                base_image.header.set_xyzt_units("mm")
                nib.save(base_image, self.base)
                image = nib.load(overlay)
                image.header.set_xyzt_units("meter")
                nib.save(image, overlay)
            with self.subTest(change=change), self.assertLogs(level="WARNING") as logs:
                self.assertEqual(self.build(nifti_in2="AnnoSplit_parental.nii.gz"), (None, 0))
            self.assertIn(expected, logs.output[0])
        self.assertFalse(list((self.project / "Report" / "Overlay").glob("*.png")))

    def test_unknown_spatial_unit_accepts_matching_geometry_without_warning_or_header_changes(self):
        overlay = self.overlay()
        for base_units, overlay_units in (("unknown", "mm"), ("mm", "unknown"),
                                          ("unknown", "unknown")):
            for path, units in ((self.base, base_units), (overlay, overlay_units)):
                image = nib.load(path)
                image.header.set_xyzt_units(units)
                nib.save(image, path)
            with self.subTest(base_units=base_units, overlay_units=overlay_units), \
                    self.assertNoLogs(level="WARNING"):
                html_path, count = self.build(nifti_in2="AnnoSplit_parental.nii.gz")
            self.assertEqual(count, 1)
            self.assertTrue(html_path.is_file())
            self.assertEqual(nib.load(self.base).header.get_xyzt_units()[0], base_units)
            self.assertEqual(nib.load(overlay).header.get_xyzt_units()[0], overlay_units)

    def test_unknown_unit_does_not_bypass_geometry_mismatch(self):
        affine = np.eye(4)
        affine[0, 3] = 0.2
        overlay = self.overlay(affine=affine)
        image = nib.load(overlay)
        image.header.set_xyzt_units("mm")
        nib.save(image, overlay)
        with self.assertLogs(level="WARNING") as logs:
            self.assertEqual(self.build(nifti_in2="AnnoSplit_parental.nii.gz"), (None, 0))
        self.assertIn("Spatial affine mismatch", logs.output[0])

    def test_small_header_differences_are_accepted_with_configurable_tolerance(self):
        affine = np.eye(4)
        affine[0, 3] = 5e-5
        affine[1, 1] += 5e-6
        self.overlay(affine=affine)
        self.assertEqual(self.build(nifti_in2="AnnoSplit_parental.nii.gz")[1], 1)
        with self.assertLogs(level="WARNING"):
            self.assertEqual(self.build(
                nifti_in2="AnnoSplit_parental.nii.gz", geometry_atol=0, geometry_rtol=0,
            ), (None, 0))

    def test_four_dimensional_input_uses_first_volume(self):
        first = np.arange(120, dtype=np.float32).reshape(4, 5, 6)
        self.save(self.base, data=np.stack([first, first + 1000], axis=-1))
        self.overlay()
        images = []
        original_imshow = Axes.imshow

        def capture(ax, data, *args, **kwargs):
            if kwargs.get("cmap") == "gray":
                images.append(data.copy())
            return original_imshow(ax, data, *args, **kwargs)

        with patch.object(Axes, "imshow", capture):
            html_path, count = self.build(nifti_in2="AnnoSplit_parental.nii.gz")
        self.assertEqual(count, 1)
        self.assertIn("(4, 5, 6, 2)", html_path.read_text())
        np.testing.assert_array_equal(images[0], np.rot90(first[:, :, 0]))

    def test_unreadable_file_is_warned_and_valid_file_still_renders(self):
        broken = self.modality / "bad.fz.fa.nii.gz"
        broken.write_bytes(b"not a nifti")
        with self.assertLogs(level="WARNING") as logs:
            html_path, count = self.build()
        self.assertEqual(count, 1)
        self.assertTrue(html_path.is_file())
        self.assertIn(broken.name, logs.output[0])

    def test_cli_creates_report_and_records_overlay_opacity(self):
        self.overlay()
        self.overlay(folder=self.modality.parent / "anat")
        csv_path = self.project / "groups.csv"
        csv_path.write_text("tDCS,Sham\nsub-NR1042,sub-NR1037\n", encoding="utf-8")
        for input1, input2, opacity, slices, modality, group in (
            ("--nifti_in1", "--nifti_in2", "--opacity", "--n-slices", "--modality", "--group"),
            ("-1", "-2", "-o", "-n", "-m", "-g"),
        ):
            with self.subTest(input1=input1), contextlib.redirect_stdout(io.StringIO()) as stdout, \
                    self.assertNoLogs(level="WARNING"):
                status = reports.main([
                    "-i", str(self.project), input1, "fz.fa.nii.gz",
                    input2, "AnnoSplit_parental.nii.gz", opacity, "60", slices, "1",
                    modality, "dwi", group, str(csv_path),
                ])
            self.assertEqual(status, 0)
            self.assertIn("(1 image(s))", stdout.getvalue())
            html = (self.project / "Report" / "Overlay" /
                    "overlay_nifti_report_fz.fa.nii.gz_AnnoSplit_parental.nii.gz.html").read_text()
            self.assertIn("Overlay opacity:</b> 60%", html)
            parsed = ReportHtml(html)
            self.assertEqual(parsed.options["groupDropdown"], [
                ("all", "All"), ("group-0", "tDCS"), ("group-1", "Sham"),
                ("unassigned", "Unassigned"),
            ])
            self.assertEqual(json.loads(parsed.entries[0]["data-groups"]), ["group-1"])

    def create_group_report(self):
        csv_path = self.project / "groups.csv"
        csv_path.write_text("tDCS,Sham\nsub-NR1042,sub-NR1037\n", encoding="utf-8")
        for subject, session, modality in (
            ("sub-NR1042", "ses-PT3", "dwi"),
            ("sub-NR1042", "ses-PT14", "anat"),
            ("sub-Unlisted", "ses-PT3", "dwi"),
        ):
            self.save(self.project / subject / session / modality / "scan.fz.fa.nii.gz")
        html_path, count = self.build(group=csv_path)
        self.assertEqual(count, 4)
        return html_path.read_text()

    def test_group_report_assigns_all_sessions_and_keeps_unlisted_subjects(self):
        parsed = ReportHtml(self.create_group_report())
        self.assertEqual(parsed.options["groupDropdown"], [
            ("all", "All"), ("group-0", "tDCS"), ("group-1", "Sham"),
            ("unassigned", "Unassigned"),
        ])
        expected = {"sub-NR1042": ["group-0"], "sub-NR1037": ["group-1"],
                    "sub-Unlisted": ["unassigned"]}
        for entry in parsed.entries:
            self.assertEqual(json.loads(entry["data-groups"]), expected[entry["data-subject"]])

    @unittest.skipUnless(shutil.which("node"), "Node.js is needed to exercise the HTML filter")
    def test_group_filter_combines_with_subject_session_modality_and_all(self):
        html = self.create_group_report()
        parsed = ReportHtml(html)
        payload = {"script": re.search(r"<script>(.*?)</script>", html, re.S).group(1),
                   "entries": parsed.entries, "cases": [
                       {"groupDropdown": "all"},
                       {"groupDropdown": "group-0"},
                       {"groupDropdown": "group-1"},
                       {"groupDropdown": "unassigned"},
                       {"groupDropdown": "group-0", "subjectDropdown": "sub-NR1042",
                        "sessionDropdown": "ses-PT14", "modalityDropdown": "anat"},
                       {"groupDropdown": "group-1", "subjectDropdown": "sub-NR1042"},
                   ]}
        javascript = """
            const fs = require('fs');
            const vm = require('vm');
            const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
            const entries = payload.entries.map(attrs => ({
                getAttribute: name => attrs[name] ?? null, style: {}, attrs,
            }));
            const controls = Object.fromEntries(
                ['groupDropdown', 'subjectDropdown', 'sessionDropdown', 'modalityDropdown']
                .map(id => [id, {value: 'all'}])
            );
            const context = {document: {
                getElementById: id => controls[id] ?? null,
                getElementsByClassName: () => entries,
                addEventListener: () => {},
            }};
            vm.createContext(context);
            vm.runInContext(payload.script, context);
            const visible = () => entries.filter(entry => entry.style.display !== 'none')
                .map(entry => entry.attrs['data-subject']).sort();
            const results = payload.cases.map(selection => {
                for (const control of Object.values(controls)) control.value = 'all';
                for (const [id, value] of Object.entries(selection)) controls[id].value = value;
                context.filterreport();
                return visible();
            });
            // Reports without CSV have no Group dropdown.
            delete controls.groupDropdown;
            for (const control of Object.values(controls)) control.value = 'all';
            context.filterreport();
            results.push(visible());
            process.stdout.write(JSON.stringify(results));
        """
        result = subprocess.run([shutil.which("node"), "-e", javascript],
                                input=json.dumps(payload), text=True, capture_output=True,
                                check=True, timeout=10)
        everyone = ["sub-NR1037", "sub-NR1042", "sub-NR1042", "sub-Unlisted"]
        self.assertEqual(json.loads(result.stdout), [
            everyone, ["sub-NR1042", "sub-NR1042"], ["sub-NR1037"], ["sub-Unlisted"],
            ["sub-NR1042"], [], everyone,
        ])

    def test_group_csv_handles_blank_cells_bom_whitespace_and_multiple_memberships(self):
        csv_path = self.project / "groups.csv"
        csv_path.write_text(
            "\ufeff tDCS , Sham \n sub-NR1042 ,sub-NR1037\nsub-NR1042,sub-NR1042\n,sub-SD624\n\n",
            encoding="utf-8",
        )
        names, memberships = reports._load_groups(csv_path)
        self.assertEqual(names, ["tDCS", "Sham"])
        self.assertEqual(memberships, {
            "sub-NR1042": ["tDCS", "Sham"], "sub-NR1037": ["Sham"], "sub-SD624": ["Sham"],
        })

    def test_group_names_are_escaped_and_do_not_conflict_with_all_filter(self):
        csv_path = self.project / "groups.csv"
        csv_path.write_text('"<Treatment & Group>",all\nsub-NR1037,sub-NR1037\n', encoding="utf-8")
        html_path, count = self.build(group=csv_path)
        self.assertEqual(count, 1)
        html = html_path.read_text()
        self.assertIn("&lt;Treatment &amp; Group&gt;", html)
        parsed = ReportHtml(html)
        self.assertEqual(parsed.options["groupDropdown"], [
            ("all", "All"), ("group-0", "<Treatment & Group>"), ("group-1", "all"),
            ("unassigned", "Unassigned"),
        ])
        self.assertEqual(json.loads(parsed.entries[0]["data-groups"]), ["group-0", "group-1"])

    def test_invalid_or_unreadable_group_csv_fails_before_creating_output(self):
        csv_path = self.project / "groups.csv"
        for text in (None, "", ",Sham\nsub-NR1042,sub-NR1037\n", "tDCS,tDCS\n",
                     "tDCS,Sham\nsub-NR1042,sub-NR1037,sub-Extra\n", 'tDCS,Sham\n"unterminated'):
            if text is not None:
                csv_path.write_text(text, encoding="utf-8")
            with self.subTest(csv=text), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    reports.main(["-i", str(self.project), "-1", "fz.fa.nii.gz", "-g", str(csv_path)])
                self.assertEqual(error.exception.code, 2)
            self.assertFalse((self.project / "Report").exists())

    def test_invalid_cli_options_fail_before_creating_output(self):
        for options in (
            ["--opacity", "60"], ["--opacity", "101"], ["--opacity", "-1"],
            ["--opacity", "nan"], ["--geometry-atol", "-1"],
            ["--geometry-rtol", "inf"], ["--n-slices", "0"],
            ["--nifti_in1", ""], ["--nifti_in2", "folder/atlas.nii.gz"],
            ["--modality", "unknown"],
        ):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    reports.main(["-i", str(self.project), "--nifti_in1", "fz.fa.nii.gz", *options])
                self.assertEqual(error.exception.code, 2)
            self.assertFalse((self.project / "Report").exists())

    def test_uncompressed_nifti_suffix_and_empty_project(self):
        self.save(self.modality / "scan_metric.nii")
        html_path, count = reports.build_display_nifti_report(
            self.project, "metric.nii", n_slices=1,
        )
        self.assertEqual(count, 1)
        self.assertEqual(html_path.name, "display_nifti_report_metric.nii.html")
        self.assertIn("scan_metric.nii", html_path.read_text())
        self.assertEqual(reports.build_display_nifti_report(self.project, "absent.nii"), (None, 0))


if __name__ == "__main__":
    unittest.main()
