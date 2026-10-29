"""Check CSV naming, ROI medians and validation before output creation."""

import contextlib
import csv
from datetime import datetime, timezone
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin/3.2.1_DTIdata_extract/DTIdata_extract.py"
spec = importlib.util.spec_from_file_location("dti_data_extract", SCRIPT)
extract = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extract)


class DTIDataExtractTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.output = self.folder / "result.csv"
        self.mapping = self.folder / "acronyms.txt"
        self.img = np.array([np.nan, 0.0001, 0.0003, 0.8]).reshape(2, 2, 1)
        self.rois = np.array([0, 1, 1, 2], dtype=np.int16).reshape(2, 2, 1)

    def write_mapping(self, text):
        self.mapping.write_text(text, encoding="utf-8")
        return self.mapping

    def read_csv(self, path=None):
        with (path or self.output).open(encoding="utf-8", newline="") as stream:
            return list(csv.reader(stream))

    def assert_invalid_preserves_output(self, img, rois, message, mapping=None, roi_file=None,
                                        stroke_mask=None):
        for existing in (False, True):
            with self.subTest(existing_output=existing):
                self.output.unlink(missing_ok=True)
                if existing:
                    self.output.write_bytes(b"previous result\n")
                with self.assertRaisesRegex(ValueError, message):
                    extract.extractDTIData(img, rois, self.output, mapping, roi_file=roi_file,
                                           stroke_mask=stroke_mask)
                if existing:
                    self.assertEqual(self.output.read_bytes(), b"previous result\n")
                else:
                    self.assertFalse(self.output.exists())

    def save(self, name, data, affine=None):
        path = self.folder / name
        nib.save(nib.Nifti1Image(data, np.eye(4) if affine is None else affine), path)
        return path

    def run_main(self, image, roi, *extra):
        process_logs_before = {path: path.read_bytes() for path in self.folder.rglob("*process.log")
                               if path.is_file()}
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = extract.main(["-i", str(image), "-r", str(roi), *map(str, extra)])
        self.assertEqual({path: path.read_bytes() for path in self.folder.rglob("*process.log")
                          if path.is_file()}, process_logs_before)
        return status, stdout.getvalue(), stderr.getvalue()

    def bids_inputs(self, metric="fa", img=None):
        dwi = self.folder / "sub-NR1037" / "ses-PT3" / "dwi"
        maps = dwi / "DSI_studio"
        maps.mkdir(parents=True, exist_ok=True)
        image = self.save(maps / f"scan.fz.{metric}.nii.gz", self.img if img is None else img)
        roi = self.save(dwi / "scan_Anno.nii.gz", self.rois)
        return dwi, image, roi

    def read_manifest(self, dwi):
        return json.loads((dwi / ".sub-NR1037_ses-PT3_dwi_aidamri_manifest.json").read_text())

    def test_output_names_preserve_image_name_and_atlas_markers(self):
        stem = "sub-NR1037_ses-PT3_dwiMR3Denoised_mcf.fz.fa"
        markers = ("Anno", "AnnoSplit", "Anno_parental", "AnnoSplit_parental")
        for extension in (".nii.gz", ".nii"):
            image = self.folder / (stem + extension)
            for marker in markers:
                for roi_marker in (marker, marker.swapcase()):
                    with self.subTest(extension=extension, roi=roi_marker):
                        roi = self.folder / f"subject_{roi_marker}_resampled{extension}"
                        actual = extract.getOutfile(roi, image)
                        self.assertEqual(Path(actual), self.folder / f"{stem}_{marker}.csv")
                        self.assertEqual(Path(extract.getOutfile(roi, image, exclude_stroke=True)),
                                         self.folder / f"{stem}_{marker}_exclude-stroke.csv")
            self.assertEqual(
                Path(extract.getOutfile(self.folder / f"custom.roi{extension}", image)),
                self.folder / f"{stem}_custom.roi.csv",
            )

    def test_csv_preserves_precision_and_quotes_labels(self):
        mapping = self.write_mapping("\n1\t\tCUL4, 5\r\n2\tfiber tracts\r\n  \n")
        img_before, rois_before = self.img.copy(), self.rois.copy()
        result = extract.extractDTIData(self.img, self.rois, self.output, mapping)
        self.assertEqual(result, self.output)
        rows = self.read_csv()
        self.assertEqual(rows[0], ["roi_id", "acronym", "median"])
        self.assertEqual(rows[1][:2], ["1", "CUL4, 5"])
        self.assertAlmostEqual(float(rows[1][2]), 0.0002)
        self.assertEqual(rows[2], ["2", "fiber tracts", "0.8"])
        self.assertEqual(len(rows), 3)
        np.testing.assert_array_equal(self.img, img_before)
        np.testing.assert_array_equal(self.rois, rois_before)

    def test_no_mapping_writes_empty_acronyms_and_skips_background(self):
        for roi_file in (None, "scan_AnnoSplit.nii.gz", "scan_AnnoSplit_parental.nii.gz"):
            with self.subTest(roi_file=roi_file):
                extract.extractDTIData(self.img, self.rois, self.output, None, roi_file=roi_file)
                rows = self.read_csv()
                self.assertEqual([row[:2] for row in rows[1:]], [["1", ""], ["2", ""]])

    def test_missing_mapping_entries_warn_and_skip_even_invalid_image_values(self):
        mapping = self.write_mapping("1\tFirst\n")
        image = self.img.copy()
        image[self.rois == 2] = np.inf
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            extract.extractDTIData(image, self.rois, self.output, mapping)
        self.assertEqual(
            stderr.getvalue(),
            f"Warning: No acronym for ROI ID 2 in '{mapping}'. Skipping region.\n",
        )
        self.assertEqual(len(self.read_csv()), 2)
        self.assertNotIn("NA", self.output.read_text())

    def test_absent_table_regions_keep_empty_rows_and_empty_atlas_without_table_has_only_header(self):
        mapping = self.write_mapping("99\tOther\n")
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            extract.extractDTIData(self.img, self.rois, self.output, mapping)
        self.assertEqual(len(stderr.getvalue().splitlines()), 2)
        expected = [["roi_id", "acronym", "median"], ["99", "Other", ""]]
        self.assertEqual(self.read_csv(), expected)
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            extract.extractDTIData(self.img, np.zeros_like(self.rois), self.output, mapping)
        self.assertEqual(stderr.getvalue(), "")
        self.assertEqual(self.read_csv(), expected)
        extract.extractDTIData(self.img, np.zeros_like(self.rois), self.output, None)
        self.assertEqual(self.read_csv(), [["roi_id", "acronym", "median"]])

    def test_table_order_and_rows_are_fixed_across_atlases_with_missing_regions(self):
        mapping = self.write_mapping("3\tThird\n2\tSecond\n0\tbackground\n1\tFirst\n2\tSecond\n")
        atlases = [self.rois, np.where(self.rois == 2, 0, self.rois), np.zeros_like(self.rois)]
        for index, rois in enumerate(atlases):
            with self.subTest(atlas=index), contextlib.redirect_stderr(io.StringIO()) as stderr:
                extract.extractDTIData(self.img, rois, self.output, mapping)
            self.assertEqual(stderr.getvalue(), "")
            rows = self.read_csv()
            self.assertEqual([row[:2] for row in rows[1:]],
                             [["3", "Third"], ["2", "Second"], ["1", "First"]])
            self.assertEqual(rows[1][2], "")
            if index == 0:
                self.assertEqual(rows[2][2], "0.8")
            else:
                self.assertEqual(rows[2][2], "")
            if index < 2:
                self.assertAlmostEqual(float(rows[3][2]), 0.0002)
            else:
                self.assertEqual(rows[3][2], "")

    def test_absent_atlas_region_is_distinct_from_region_fully_excluded_by_stroke(self):
        mapping = self.write_mapping("3\tMissing\n2\tSecond\n1\tFirst\n")
        mask = np.array([0, 0, 0, 1]).reshape(self.img.shape)
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            extract.extractDTIData(self.img, self.rois, self.output, mapping, stroke_mask=mask)
        rows = self.read_csv()
        self.assertEqual(rows[1], ["3", "Missing", ""])
        self.assertEqual(rows[2], ["2", "Second", ""])
        self.assertAlmostEqual(float(rows[3][2]), 0.0002)
        self.assertEqual(len(stderr.getvalue().splitlines()), 1)
        self.assertIn("All voxels in ROI ID 2 excluded by stroke mask", stderr.getvalue())
        self.assertNotIn("ROI ID 3", stderr.getvalue())

    def test_shapes_must_be_matching_and_three_dimensional(self):
        cases = [
            (self.img[..., None], self.rois, "must be 3D"),
            (self.img, self.rois[..., None], "must be 3D"),
            (self.img[..., None], self.rois[..., None], "must be 3D"),
            (self.img, np.ones((3, 2, 1)), "do not match"),
        ]
        for img, rois, message in cases:
            with self.subTest(image_shape=img.shape, roi_shape=rois.shape):
                self.assert_invalid_preserves_output(img, rois, message)

    def test_roi_labels_must_be_finite_integers(self):
        for value in (np.nan, np.inf, -np.inf, 1.5):
            rois = self.rois.astype(float)
            rois[0, 0, 0] = value
            with self.subTest(label=value):
                self.assert_invalid_preserves_output(self.img, rois, "ROI labels must")
        extract.extractDTIData(self.img, self.rois.astype(float), self.output, None)
        self.assertEqual(len(self.read_csv()), 3)

    def test_nonfinite_values_in_output_regions_preserve_existing_output(self):
        for value in (np.nan, np.inf, -np.inf):
            image = self.img.copy()
            # The error is in the second ROI, after the first median is collected.
            image[self.rois == 2] = value
            with self.subTest(value=value):
                self.assert_invalid_preserves_output(image, self.rois, "ROI ID 2")

    def test_stroke_mask_excludes_voxels_without_changing_inputs(self):
        mask = np.array([1, 0, 1, 0], dtype=np.uint8).reshape(self.img.shape)
        originals = [data.copy() for data in (self.img, self.rois, mask)]
        extract.extractDTIData(self.img, self.rois, self.output, None, stroke_mask=mask)
        rows = self.read_csv()
        self.assertAlmostEqual(float(rows[1][2]), 0.0001)
        self.assertEqual(rows[2], ["2", "", "0.8"])
        for data, before in zip((self.img, self.rois, mask), originals):
            np.testing.assert_array_equal(data, before)

    def test_zero_stroke_mask_preserves_medians(self):
        extract.extractDTIData(self.img, self.rois, self.output, None)
        expected = self.read_csv()
        extract.extractDTIData(self.img, self.rois, self.output, None,
                               stroke_mask=np.zeros_like(self.img))
        self.assertEqual(self.read_csv(), expected)

    def test_fully_excluded_rois_keep_ids_acronyms_and_empty_medians(self):
        mapping = self.write_mapping("1\tFirst\n2\tSecond\n")
        for mask in ((0, 0, 0, 1), (1, 1, 1, 1)):
            with self.subTest(mask=mask), contextlib.redirect_stderr(io.StringIO()) as stderr:
                extract.extractDTIData(self.img, self.rois, self.output, mapping,
                                       stroke_mask=np.array(mask).reshape(self.img.shape))
            rows = self.read_csv()
            self.assertEqual(rows[1][:2], ["1", "First"])
            self.assertEqual(rows[2], ["2", "Second", ""])
            self.assertIn("All voxels in ROI ID 2 excluded by stroke mask", stderr.getvalue())
            self.assertIn("[subject=unknown, session=unknown, metric=unknown]", stderr.getvalue())
            if all(mask):
                self.assertEqual(rows[1], ["1", "First", ""])
                self.assertIn("All voxels in ROI ID 1 excluded by stroke mask", stderr.getvalue())
            else:
                self.assertAlmostEqual(float(rows[1][2]), 0.0002)
            self.assertNotIn("nan", self.output.read_text().lower())

    def test_invalid_stroke_masks_preserve_output(self):
        cases = [
            (np.zeros((*self.img.shape, 1)), "Stroke mask must be 3D"),
            (np.zeros((3, 2, 1)), "stroke mask shape.*do not match"),
        ]
        for value in (-1, 2, 255, 0.5, np.nan, np.inf, -np.inf):
            mask = np.zeros_like(self.img)
            mask[0, 0, 0] = value
            cases.append((mask, "Stroke mask must be binary"))
        for mask, message in cases:
            with self.subTest(mask=mask):
                self.assert_invalid_preserves_output(self.img, self.rois, message,
                                                     stroke_mask=mask)

    def test_nonfinite_metric_values_are_allowed_only_in_excluded_voxels(self):
        mask = np.array([0, 0, 1, 0]).reshape(self.img.shape)
        for value in (np.nan, np.inf, -np.inf):
            image = self.img.copy()
            image[1, 0, 0] = value
            with self.subTest(value=value):
                extract.extractDTIData(image, self.rois, self.output, None, stroke_mask=mask)
                self.assertAlmostEqual(float(self.read_csv()[1][2]), 0.0001)
                image[0, 1, 0] = value
                self.assert_invalid_preserves_output(image, self.rois, "ROI ID 1",
                                                     stroke_mask=mask)

    def test_cli_stroke_mask_suffix_log_and_manifest_preserve_inputs_and_regular_csv(self):
        dwi, image, roi = self.bids_inputs()
        mask = self.save(dwi / "stroke mask.nii.gz",
                         np.array([0, 0, 1, 1], dtype=np.uint8).reshape(self.img.shape))
        original_inputs = {path: path.read_bytes() for path in (image, roi, mask)}
        regular_output = Path(extract.getOutfile(roi, image))
        regular_output.write_bytes(b"regular export\n")
        output = image.parent / "scan.fz.fa_Anno_exclude-stroke.csv"
        for _ in range(2):
            # Exercise relative paths, including spaces in the mask filename.
            relative_mask = os.path.relpath(mask)
            status, stdout, stderr = self.run_main(image, roi, "--exclude-stroke-mask", relative_mask)
            self.assertEqual(status, 0, stderr)
            self.assertIn(f"Done: {output}", stdout)
            self.assertIn(f"Exclude stroke mask: {mask}", stdout)
            self.assertIn("All voxels in ROI ID 2 excluded by stroke mask", stderr)
            self.assertIn("Warning: [subject=sub-NR1037, session=ses-PT3, metric=FA] "
                          "All voxels in ROI ID 2 excluded by stroke mask. Leaving median empty.", stderr)
            rows = self.read_csv(output)
            self.assertAlmostEqual(float(rows[1][2]), 0.0001)
            self.assertEqual(rows[2], ["2", "", ""])
        log = extract._analysis_log_path(image).read_text()
        self.assertIn(f"Exclude stroke mask: {mask}", log)
        self.assertIn("All voxels in ROI ID 2 excluded by stroke mask", log)
        self.assertIn("Warning: [subject=sub-NR1037, session=ses-PT3, metric=FA] "
                      "All voxels in ROI ID 2 excluded by stroke mask. Leaving median empty.", log)
        self.assertEqual(regular_output.read_bytes(), b"regular export\n")
        stage = self.read_manifest(dwi)["stages"]["extraction"]
        self.assertIn(output.relative_to(dwi).as_posix(), stage["created_files"])
        self.assertIn(output.relative_to(dwi).as_posix(), stage["modified_files"])
        for path, before in original_inputs.items():
            self.assertEqual(path.read_bytes(), before)

    def test_cli_zero_stroke_mask_still_uses_suffix(self):
        image = self.save("scan.fa.nii", self.img)
        roi = self.save("scan_Anno.nii", self.rois)
        mask = self.save("stroke.nii", np.zeros_like(self.rois))
        status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        regular = Path(extract.getOutfile(roi, image))
        status, _, stderr = self.run_main(image, roi, "--exclude-stroke-mask", mask)
        self.assertEqual(status, 0, stderr)
        self.assertEqual((self.folder / "scan.fa_Anno_exclude-stroke.csv").read_bytes(), regular.read_bytes())

    def test_cli_stroke_mask_validation_preserves_outputs(self):
        image = self.save("scan.fa.nii.gz", self.img)
        roi = self.save("scan_Anno.nii.gz", self.rois)
        shifted = np.eye(4)
        shifted[0, 3] = 0.00002
        cases = [
            (self.folder / "missing.nii.gz", "No such file"),
            (self.save("wrong_shape.nii.gz", np.zeros((3, 2, 1))), "do not match"),
            (self.save("4d.nii.gz", np.zeros((*self.img.shape, 1))), "must be 3D"),
            (self.save("invalid.nii.gz", np.full_like(self.img, np.nan)), "must be binary"),
            (self.save("shifted.nii.gz", np.zeros_like(self.img), shifted), "affines do not match"),
        ]
        regular = Path(extract.getOutfile(roi, image))
        regular.write_bytes(b"regular export\n")
        output = self.folder / "scan.fa_Anno_exclude-stroke.csv"
        for mask, message in cases:
            for existing in (False, True):
                with self.subTest(mask=mask.name, existing=existing):
                    output.unlink(missing_ok=True)
                    if existing:
                        output.write_bytes(b"previous stroke export\n")
                    status, stdout, stderr = self.run_main(image, roi, "--exclude-stroke-mask", mask)
                    self.assertEqual(status, 1)
                    self.assertIn(message, stderr)
                    self.assertNotIn("Done:", stdout)
                    self.assertIn(message, extract._analysis_log_path(image).read_text())
                    if existing:
                        self.assertEqual(output.read_bytes(), b"previous stroke export\n")
                    else:
                        self.assertFalse(output.exists())
                    self.assertEqual(regular.read_bytes(), b"regular export\n")

    def test_cli_accepts_stroke_mask_affine_within_tolerance(self):
        image = self.save("scan.fa.nii.gz", self.img)
        roi = self.save("scan_Anno.nii.gz", self.rois)
        affine = np.eye(4)
        affine[0, 3] = 0.000005
        mask = self.save("stroke.nii.gz", np.zeros_like(self.img), affine)
        status, _, stderr = self.run_main(image, roi, "--exclude-stroke-mask", mask)
        self.assertEqual(status, 0, stderr)

    def test_bad_translation_lines_report_location_and_preserve_output(self):
        for line, message in (
            ("1", "Missing ROI ID or label"),
            ("invalid\tLabel", "Invalid integer ROI ID"),
            ("1\tDifferent", "Conflicting labels for ROI ID 1"),
        ):
            mapping = self.write_mapping("1\tFirst\n" + line + "\n")
            with self.subTest(line=line):
                self.assert_invalid_preserves_output(self.img, self.rois, message, mapping)
                with self.assertRaises(ValueError) as error:
                    extract._read_mapping(mapping)
                self.assertIn(f"'{mapping}', line 2", str(error.exception))
        self.write_mapping("1\tFirst\n1\tFirst\n")
        self.assertEqual(extract._read_mapping(self.mapping), {1: "First"})

    def test_both_repository_translation_files_are_supported(self):
        ara = extract._read_mapping(ROOT / "lib/acronyms_ARA.txt")
        fmri = extract._read_mapping(ROOT / "lib/acronym_rsfMRI.txt")
        self.assertEqual(ara[1009], "fiber tracts")
        self.assertEqual(fmri[22], "PTLp")

    def test_split_atlases_write_distinct_hemisphere_labels_and_medians(self):
        mapping = self.write_mapping("22\tPTLp\n")
        before = mapping.read_bytes()
        rois = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        for marker in ("AnnoSplit", "AnnoSplit_parental", "aNnOsPlIt_PaReNtAl"):
            with self.subTest(marker=marker), contextlib.redirect_stderr(io.StringIO()) as stderr:
                extract.extractDTIData(
                    self.img, rois, self.output, mapping, roi_file=f"subject_{marker}.nii.gz"
                )
                rows = self.read_csv()
                self.assertEqual([row[:2] for row in rows[1:]], [["22", "L_PTLp"], ["2022", "R_PTLp"]])
                self.assertAlmostEqual(float(rows[1][2]), 0.0002)
                self.assertAlmostEqual(float(rows[2][2]), 0.8)
                self.assertEqual(stderr.getvalue(), "")
        self.assertEqual(mapping.read_bytes(), before)

    def test_unsplit_roi_names_leave_mapping_unchanged(self):
        mapping = self.write_mapping("1\tFirst\n2\tSecond\n")
        for roi_file in (None, "subject_Anno.nii.gz", "subject_Anno_parental.nii.gz"):
            with self.subTest(roi_file=roi_file):
                extract.extractDTIData(self.img, self.rois, self.output, mapping, roi_file=roi_file)
                self.assertEqual(
                    [row[:2] for row in self.read_csv()[1:]], [["1", "First"], ["2", "Second"]]
                )

    def test_split_order_is_left_then_right_and_prelabelled_order_is_preserved(self):
        rois = np.array([22, 31, 2031, 2022], dtype=np.int16).reshape(self.img.shape)
        image = np.array([0.7, 0.1, 0.2, 0.9]).reshape(self.img.shape)
        cases = [
            ("31\tACA\n22\tPTLp\n",
             [["31", "L_ACA", "0.1"], ["22", "L_PTLp", "0.7"],
              ["2031", "R_ACA", "0.2"], ["2022", "R_PTLp", "0.9"]]),
            ("2022\tR_PTLp\n31\tL_ACA\n2031\tR_ACA\n22\tL_PTLp\n",
             [["2022", "R_PTLp", "0.9"], ["31", "L_ACA", "0.1"],
              ["2031", "R_ACA", "0.2"], ["22", "L_PTLp", "0.7"]]),
        ]
        for text, expected in cases:
            for atlas in ("scan_AnnoSplit.nii.gz", "scan_AnnoSplit_parental.nii.gz"):
                with self.subTest(text=text, atlas=atlas):
                    mapping = self.write_mapping(text)
                    extract.extractDTIData(image, rois, self.output, mapping, roi_file=atlas)
                    self.assertEqual(self.read_csv()[1:], expected)

    def test_large_base_ids_are_expanded_without_hemisphere_threshold(self):
        base_id = 304325711
        mapping = self.write_mapping(f"{base_id}\tretina\n")
        rois = np.array([0, base_id, base_id, base_id + 2000], dtype=np.int32).reshape(2, 2, 1)
        extract.extractDTIData(self.img, rois, self.output, mapping, roi_file="scan_AnnoSplit.nii.gz")
        self.assertEqual(
            [row[:2] for row in self.read_csv()[1:]],
            [[str(base_id), "L_retina"], [str(base_id + 2000), "R_retina"]],
        )

    def test_split_mapping_does_not_expand_background(self):
        mapping = self.write_mapping("0\tbackground\n22\tPTLp\n")
        self.assertEqual(extract._read_mapping(mapping, True), {22: "L_PTLp", 2022: "R_PTLp"})
        self.write_mapping("0\tbackground\n22\tL_PTLp\n2022\tR_PTLp\n")
        self.assertEqual(extract._read_mapping(mapping, True), {22: "L_PTLp", 2022: "R_PTLp"})

    def test_repository_mappings_expand_and_existing_split_table_is_preserved(self):
        for name in ("acronyms_ARA.txt", "acronym_rsfMRI.txt"):
            with self.subTest(name=name):
                source = ROOT / "lib" / name
                base = extract._read_mapping(source)
                expanded = extract._read_mapping(source, True)
                self.assertEqual(len(expanded), 2 * len(base))
                for roi_id, acronym in base.items():
                    self.assertEqual(expanded[roi_id], "L_" + acronym)
                    self.assertEqual(expanded[roi_id + 2000], "R_" + acronym)
        split_file = SCRIPT.parent / "acronyms_splitted_ARA.txt"
        original = extract._read_mapping(split_file)
        self.assertEqual(extract._read_mapping(split_file, True), original)
        self.assertEqual(original, extract._read_mapping(ROOT / "lib/acronyms_ARA.txt", True))
        self.assertEqual(list(original.items()),
                         list(extract._read_mapping(ROOT / "lib/acronyms_ARA.txt", True).items()))

    def test_existing_split_mapping_produces_no_double_prefixes(self):
        mapping = self.write_mapping("22\tL_PTLp\n2022\tR_PTLp\n")
        rois = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        extract.extractDTIData(self.img, rois, self.output, mapping, roi_file="scan_AnnoSplit.nii.gz")
        self.assertEqual(
            [row[:2] for row in self.read_csv()[1:]], [["22", "L_PTLp"], ["2022", "R_PTLp"]]
        )

    def test_conflicting_split_mappings_preserve_output(self):
        for text, message in (
            ("22\tL_PTLp\n31\tACA\n", "Mixed hemisphere-labelled"),
            ("22\tPTLp\n2022\tOther\n", "Hemisphere ID collision for ROI ID 2022"),
            ("2022\tOther\n22\tPTLp\n", "Hemisphere ID collision for ROI ID 2022"),
        ):
            with self.subTest(mapping=text):
                mapping = self.write_mapping(text)
                self.assert_invalid_preserves_output(
                    self.img, self.rois, message, mapping, roi_file="subject_AnnoSplit_parental.nii.gz"
                )

    def test_unknown_ids_in_split_atlas_warn_and_skip(self):
        mapping = self.write_mapping("22\tPTLp\n")
        rois = np.array([0, 22, 2022, 9999], dtype=np.int16).reshape(2, 2, 1)
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            extract.extractDTIData(self.img, rois, self.output, mapping, roi_file="scan_AnnoSplit.nii.gz")
        self.assertIn("No acronym for ROI ID 9999", stderr.getvalue())
        self.assertEqual(len(stderr.getvalue().splitlines()), 1)
        self.assertEqual(len(self.read_csv()), 3)

    def test_cli_exports_nifti_inputs_with_mapping(self):
        image = self.save("subject.fz.fa.nii.gz", self.img)
        mapping = self.write_mapping("22\tPTLp\n")
        rois = np.array([0, 22, 22, 2022], dtype=np.int16).reshape(2, 2, 1)
        for marker in ("AnnoSplit", "AnnoSplit_parental"):
            with self.subTest(marker=marker):
                roi = self.save(f"subject_{marker}.nii.gz", rois)
                status, stdout, stderr = self.run_main(image, roi, "-t", mapping)
                output = self.folder / f"subject.fz.fa_{marker}.csv"
                self.assertEqual(status, 0, stderr)
                self.assertIn(f"Done: {output}", stdout)
                self.assertEqual(stderr.count("Manifest tracking skipped"), 1)
                self.assertIn(f"Done: {output}", (self.folder / "analysis.log").read_text())
                self.assertEqual(
                    [row[:2] for row in self.read_csv(output)[1:]], [["22", "L_PTLp"], ["2022", "R_PTLp"]]
                )

    def test_cli_requires_image_and_roi_options(self):
        for args in ([], ["-i", "image.nii.gz"], ["-r", "roi.nii.gz"],
                     ["-i", "image.nii.gz", "-r", "roi.nii.gz", "--exclude-stroke-mask"]):
            with self.subTest(args=args), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    extract.main(args)
                self.assertEqual(error.exception.code, 2)

    def test_cli_affine_tolerance_and_failure_preserve_output(self):
        image = self.save("subject.fz.fa.nii", self.img)
        for shift, succeeds in ((0.000005, True), (0.00002, False)):
            affine = np.eye(4)
            affine[0, 3] = shift
            roi = self.save("subject_Anno.nii", self.rois, affine)
            output = Path(extract.getOutfile(roi, image))
            for existing in (False, True):
                with self.subTest(shift=shift, existing=existing):
                    output.unlink(missing_ok=True)
                    if existing:
                        output.write_bytes(b"previous result\n")
                    status, _, stderr = self.run_main(image, roi)
                    if succeeds:
                        self.assertEqual(status, 0, stderr)
                        self.assertEqual(len(self.read_csv(output)), 3)
                    else:
                        self.assertEqual(status, 1)
                        self.assertIn("affines do not match", stderr)
                        self.assertIn("affines do not match", (self.folder / "analysis.log").read_text())
                        if existing:
                            self.assertEqual(output.read_bytes(), b"previous result\n")
                        else:
                            self.assertFalse(output.exists())

    def test_command_line_exit_status_and_paths_with_spaces(self):
        image = self.save("image with spaces.fz.fa.nii.gz", self.img)
        roi = self.save("roi Anno.nii.gz", self.rois)
        command = [sys.executable, "-B", str(SCRIPT), "--image_file", str(image), "--roi_file", str(roi)]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Done:", result.stdout)
        output = Path(extract.getOutfile(roi, image))
        previous = output.read_bytes()
        image_data = self.img.copy()
        image_data[self.rois == 1] = np.nan
        self.save(image.name, image_data)
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 1)
        self.assertNotIn("Done:", result.stdout)
        self.assertEqual(result.stderr.count("ROI ID 1"), 1)
        log = (self.folder / "analysis.log").read_text()
        self.assertIn("ROI ID 1", log)
        self.assertNotIn("Traceback", log)
        self.assertEqual(output.read_bytes(), previous)
        self.assertFalse((self.folder / "process.log").exists())

    def test_cli_appends_success_warnings_and_errors_to_bids_analysis_log(self):
        dwi = self.folder / "sub-NR1037" / "ses-PT3" / "dwi"
        maps = dwi / "DSI_studio"
        maps.mkdir(parents=True)
        image = self.save(maps / "scan.fz.fa.nii.gz", self.img)
        labels = np.array([0, 22, 2022, 9999], dtype=np.int16).reshape(2, 2, 1)
        roi = self.save(dwi / "scan_AnnoSplit_parental.nii.gz", labels)
        translator = self.write_mapping("22\tPTLp\n")
        logfile = dwi / "sub-NR1037_ses-PT3_dwi_analysis.log"
        logfile.write_text("Existing analysis entries\n", encoding="utf-8")
        process_log = dwi / "sub-NR1037_ses-PT3_dwi_process.log"
        process_log.write_bytes(b"Existing pipeline entries\n")
        git_information = "AIDAmri test Git information\nCommit: test-commit\n"
        (self.folder / "AIDAmri_git_information.txt").write_text(git_information, encoding="utf-8")
        build_path = patch("common.script_logging.BUILD_GIT_INFO_PATH", self.folder / "missing-build-info")
        build_path.start()
        self.addCleanup(build_path.stop)
        original_stdout, original_stderr = sys.stdout, sys.stderr
        for _ in range(2):
            status, stdout, stderr = self.run_main(image, roi, "-t", translator)
            self.assertEqual(status, 0)
            self.assertEqual(stdout.count("Done:"), 1)
            self.assertNotIn("Warning:", stdout)
            self.assertEqual(stderr.count("No acronym for ROI ID 9999"), 1)
            self.assertNotIn("Done:", stderr)
            self.assertIs(sys.stdout, original_stdout)
            self.assertIs(sys.stderr, original_stderr)
        output = Path(extract.getOutfile(roi, image))
        previous = output.read_bytes()
        invalid = self.img.copy()
        invalid[labels == 22] = np.nan
        self.save(image, invalid)
        status, stdout, stderr = self.run_main(image, roi, "-t", translator)
        self.assertEqual(status, 1)
        self.assertNotIn("Done:", stdout)
        self.assertEqual(stderr.count("Image contains NaN or infinity in ROI ID 22"), 1)
        self.assertIs(sys.stdout, original_stdout)
        self.assertIs(sys.stderr, original_stderr)
        self.assertEqual(output.read_bytes(), previous)
        log = logfile.read_text(encoding="utf-8")
        self.assertTrue(log.startswith("Existing analysis entries\n"))
        self.assertNotIn("Existing pipeline entries", log)
        self.assertEqual(process_log.read_bytes(), b"Existing pipeline entries\n")
        self.assertEqual(log.count("=== DTIdata_extract.py "), 3)
        self.assertRegex(log, r"DTIdata_extract\.py \d{2} [A-Za-z]+ \d{4} \d{2}:\d{2}:\d{2} CE(?:S)?T ===")
        self.assertEqual(log.count("Done:"), 2)
        # Unmapped atlas IDs are reported before processing the table, including failed runs.
        self.assertEqual(log.count("No acronym for ROI ID 9999"), 3)
        self.assertIn("Image contains NaN or infinity in ROI ID 22", log)
        self.assertIn(f"Image: {image}", log)
        self.assertIn(f"ROI: {roi}", log)
        self.assertIn(f"Translation: {translator}", log)
        self.assertEqual(list(maps.glob("*analysis.log")), [])
        self.assertEqual(log.count(git_information), 3)
        self.assertTrue(log.endswith(git_information))
        self.assertNotIn(git_information, stdout + stderr)

    def test_cli_uses_berlin_timestamps_in_summer_and_winter(self):
        _, image, roi = self.bids_inputs()
        cases = (
            (datetime(2026, 10, 16, 13, 15, 41, tzinfo=timezone.utc), "16 October 2026 15:15:41 CEST"),
            (datetime(2026, 1, 16, 23, 15, 41, tzinfo=timezone.utc), "17 January 2026 00:15:41 CET"),
        )
        for instant, expected in cases:
            with self.subTest(instant=instant), patch("common.script_logging.datetime") as clock:
                clock.now.side_effect = lambda tz=None: instant.astimezone(tz)
                status, stdout, stderr = self.run_main(image, roi)
            self.assertEqual(status, 0, stderr)
            header = f"=== DTIdata_extract.py {expected} ==="
            self.assertIn(header, stdout.splitlines())
            self.assertIn(header, extract._analysis_log_path(image).read_text().splitlines())

    def test_only_matching_parent_header_marker_suppresses_the_header(self):
        dwi, image, roi = self.bids_inputs()
        for marker, suppressed in ((str(dwi), True), (str(self.folder), False), (".", False), ("", False)):
            before = extract._analysis_log_path(image).read_text() if extract._analysis_log_path(image).exists() else ""
            with self.subTest(marker=marker), patch.dict(os.environ, {extract.PARENT_HEADER_LOG_ENV: marker}):
                status, stdout, stderr = self.run_main(image, roi)
            self.assertEqual(status, 0, stderr)
            self.assertEqual(stdout.count("=== DTIdata_extract.py "), 0 if suppressed else 1)
            added = extract._analysis_log_path(image).read_text()[len(before):]
            self.assertEqual(added.count("=== DTIdata_extract.py "), 0 if suppressed else 1)
            self.assertIn("Image:", added)
            self.assertIn("Done:", added)
            # Header ownership must not suppress standalone manifest tracking.
            self.assertTrue(self.read_manifest(dwi)["stages"]["extraction"]["created_files"])

    def test_git_block_follows_manifest_errors(self):
        dwi, image, roi = self.bids_inputs()
        git_information = "Git information after tracking failure\n"
        source = self.folder / "AIDAmri_git_information.txt"
        source.write_text(git_information)
        with patch("common.script_logging.BUILD_GIT_INFO_PATH", source), \
                patch("common.artifact_manifest._merge_changes", side_effect=OSError("manifest denied")):
            status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 1)
        self.assertIn("manifest denied", stderr)
        log = extract._analysis_log_path(image).read_text()
        self.assertEqual(log.count(git_information), 1)
        self.assertTrue(log.endswith(git_information))

    def test_git_logging_has_a_separate_folder_scoped_parent_marker(self):
        dwi, image, roi = self.bids_inputs()
        # Owning Git logging must not suppress the standalone manifest tracker.
        for marker, suppressed in ((str(dwi), True), (str(self.folder), False), (".", False), ("", False)):
            with self.subTest(marker=marker), patch.dict(os.environ, {extract.PARENT_GIT_LOG_ENV: marker}), \
                    patch.object(extract, "append_git_information") as append:
                status, _, stderr = self.run_main(image, roi)
            self.assertEqual(status, 0, stderr)
            self.assertEqual(append.call_count, 0 if suppressed else 1)
            self.assertTrue(self.read_manifest(dwi)["stages"]["extraction"]["created_files"])

    def test_git_write_failure_is_reported_and_keeps_csv(self):
        dwi, image, roi = self.bids_inputs()
        with patch.object(extract, "append_git_information", side_effect=OSError("Git log denied")):
            status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 1)
        self.assertIn("Cannot append Git information", stderr)
        self.assertIn("Git log denied", stderr)
        self.assertIn(stderr.strip(), extract._analysis_log_path(image).read_text())
        self.assertTrue(Path(extract.getOutfile(roi, image)).is_file())

    def test_missing_git_information_writes_warning_to_log_only(self):
        dwi, image, roi = self.bids_inputs()
        with patch("common.script_logging.BUILD_GIT_INFO_PATH", self.folder / "absent"):
            status, stdout, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        self.assertNotIn("Could not read AIDAmri Git information", stdout + stderr)
        log = extract._analysis_log_path(image).read_text()
        self.assertEqual(log.count("WARNING: Could not read AIDAmri Git information"), 1)
        self.assertIn("WARNING: Could not read AIDAmri Git information", log.splitlines()[-1])

    def test_unwritable_log_aborts_before_csv_is_opened(self):
        image = self.save("scan.fa.nii.gz", self.img)
        roi = self.save("scan_Anno.nii.gz", self.rois)
        (self.folder / "analysis.log").mkdir()
        output = Path(extract.getOutfile(roi, image))
        output.write_bytes(b"previous result\n")
        status, stdout, stderr = self.run_main(image, roi)
        self.assertEqual(status, 1)
        self.assertEqual(stdout, "")
        self.assertIn("Cannot write analysis log", stderr)
        self.assertEqual(output.read_bytes(), b"previous result\n")

    def test_new_csv_and_log_are_tracked_under_extraction_not_inputs(self):
        dwi, image, roi = self.bids_inputs()
        original_image, original_roi = image.read_bytes(), roi.read_bytes()
        status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        self.assertEqual(stderr, "")
        manifest = self.read_manifest(dwi)
        self.assertEqual(manifest["mode"], "dwi")
        self.assertEqual(manifest["tracking_started_before_stage"], "extraction")
        self.assertEqual(set(manifest["stages"]), {"extraction"})
        stage = manifest["stages"]["extraction"]
        self.assertEqual(stage["created_files"], [
            "DSI_studio/scan.fz.fa_Anno.csv", "sub-NR1037_ses-PT3_dwi_analysis.log",
        ])
        self.assertEqual(stage["modified_files"], [])
        self.assertEqual(stage["created_directories"], [])
        self.assertEqual(image.read_bytes(), original_image)
        self.assertEqual(roi.read_bytes(), original_roi)
        self.assertEqual(list(image.parent.glob(".*aidamri_manifest.json")), [])

    def test_existing_outputs_are_modified_and_other_stages_preserved(self):
        dwi, image, roi = self.bids_inputs()
        registration = extract.start_output_tracking(dwi, "dwi", "registration")
        (dwi / "registration_result.txt").write_text("registered")
        registration.finalize()
        previous_registration = self.read_manifest(dwi)["stages"]["registration"]
        output = Path(extract.getOutfile(roi, image))
        output.write_text("preexisting CSV")
        logfile = dwi / "sub-NR1037_ses-PT3_dwi_analysis.log"
        logfile.write_text("existing pipeline log\n")
        status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        manifest = self.read_manifest(dwi)
        self.assertEqual(manifest["stages"]["registration"], previous_registration)
        stage = manifest["stages"]["extraction"]
        self.assertEqual(stage["created_files"], [])
        self.assertEqual(stage["modified_files"], [
            "DSI_studio/scan.fz.fa_Anno.csv", logfile.name,
        ])
        self.assertTrue(logfile.read_text().startswith("existing pipeline log\n"))

    def test_legacy_processing_csv_is_not_migrated_to_extraction_ownership(self):
        dwi, image, roi = self.bids_inputs()
        legacy_tracker = extract.start_output_tracking(dwi, "dwi", "processing")
        output = Path(extract.getOutfile(roi, image))
        output.write_text("previous export")
        legacy_tracker.finalize()
        previous = self.read_manifest(dwi)["stages"]["processing"]
        status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        manifest = self.read_manifest(dwi)
        self.assertEqual(manifest["stages"]["processing"], previous)
        relative = output.relative_to(dwi).as_posix()
        self.assertIn(relative, manifest["stages"]["processing"]["created_files"])
        self.assertNotIn(relative, manifest["stages"]["extraction"]["created_files"])
        self.assertIn(relative, manifest["stages"]["extraction"]["modified_files"])

    def test_old_processing_parent_marker_does_not_suppress_extraction_tracking(self):
        dwi, image, roi = self.bids_inputs()
        with patch.dict(os.environ, {"AIDAMRI_DTI_PROCESSING_TRACKED_BY_PARENT": str(dwi)}), \
                patch.object(extract, "start_output_tracking", wraps=extract.start_output_tracking) as start:
            status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        start.assert_called_once_with(dwi, "dwi", "extraction")
        self.assertEqual(set(self.read_manifest(dwi)["stages"]), {"extraction"})

    def test_multiple_metrics_and_repeated_main_calls_merge_manifest_immediately(self):
        dwi, image, roi = self.bids_inputs()
        status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        self.assertIn("DSI_studio/scan.fz.fa_Anno.csv", self.read_manifest(dwi)["stages"]["extraction"]["created_files"])
        _, ad_image, roi = self.bids_inputs("ad")
        status, _, stderr = self.run_main(ad_image, roi)
        self.assertEqual(status, 0, stderr)
        values = self.img.copy()
        values[self.rois == 1] = 0.6
        self.save(image, values)
        status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        stage = self.read_manifest(dwi)["stages"]["extraction"]
        self.assertEqual(stage["created_files"], [
            "DSI_studio/scan.fz.ad_Anno.csv", "DSI_studio/scan.fz.fa_Anno.csv",
            "sub-NR1037_ses-PT3_dwi_analysis.log",
        ])
        self.assertIn("DSI_studio/scan.fz.fa_Anno.csv", stage["modified_files"])
        self.assertIn("sub-NR1037_ses-PT3_dwi_analysis.log", stage["modified_files"])

    def test_failed_extraction_tracks_log_without_claiming_csv(self):
        invalid = self.img.copy()
        invalid[self.rois == 1] = np.nan
        dwi, image, roi = self.bids_inputs(img=invalid)
        output = Path(extract.getOutfile(roi, image))
        for existing in (False, True):
            with self.subTest(existing=existing):
                if existing:
                    output.write_text("previous result\n")
                status, _, stderr = self.run_main(image, roi)
                self.assertEqual(status, 1)
                self.assertIn("ROI ID 1", stderr)
                stage = self.read_manifest(dwi)["stages"]["extraction"]
                self.assertNotIn("DSI_studio/scan.fz.fa_Anno.csv", stage["created_files"])
                self.assertNotIn("DSI_studio/scan.fz.fa_Anno.csv", stage["modified_files"])
                if existing:
                    self.assertEqual(output.read_text(), "previous result\n")
                else:
                    self.assertFalse(output.exists())

    def test_tracking_start_failure_is_reported_in_terminal_and_log(self):
        dwi, image, roi = self.bids_inputs()
        with patch.object(extract, "start_output_tracking", side_effect=PermissionError("tracking denied")) as start:
            status, stdout, stderr = self.run_main(image, roi)
        start.assert_called_once_with(dwi, "dwi", "extraction")
        self.assertEqual(status, 1)
        self.assertNotIn("Done:", stdout)
        self.assertIn("Cannot start extraction manifest tracking", stderr)
        self.assertIn("tracking denied", stderr)
        self.assertFalse(Path(extract.getOutfile(roi, image)).exists())
        self.assertIn(stderr.strip(), (dwi / "sub-NR1037_ses-PT3_dwi_analysis.log").read_text())

    def test_manifest_write_failure_keeps_csv_and_returns_failure(self):
        dwi, image, roi = self.bids_inputs()
        with patch("common.artifact_manifest._merge_changes", side_effect=OSError("manifest denied")):
            status, _, stderr = self.run_main(image, roi)
        self.assertEqual(status, 1)
        self.assertIn("Cannot finalize extraction manifest", stderr)
        self.assertIn("manifest denied", stderr)
        self.assertTrue(Path(extract.getOutfile(roi, image)).is_file())
        log = (dwi / "sub-NR1037_ses-PT3_dwi_analysis.log").read_text()
        self.assertIn(stderr.strip(), log)

    def test_non_bids_inputs_warn_and_continue_without_manifest(self):
        for maps in (self.folder, self.folder / "animal" / "time" / "dwi" / "DSI_studio"):
            with self.subTest(maps=maps):
                maps.mkdir(parents=True, exist_ok=True)
                image = self.save(maps / "scan.fa.nii.gz", self.img)
                roi = self.save(maps / "scan_Anno.nii.gz", self.rois)
                status, _, stderr = self.run_main(image, roi)
                self.assertEqual(status, 0, stderr)
                self.assertEqual(stderr.count("Manifest tracking skipped"), 1)
                log = extract._analysis_log_path(image).read_text()
                self.assertIn(stderr.strip(), log)
                self.assertTrue(Path(extract.getOutfile(roi, image)).is_file())
                self.assertEqual(list(self.folder.rglob("*aidamri_manifest.json")), [])

    def test_matching_parent_folder_disables_only_manifest_tracking(self):
        dwi, image, roi = self.bids_inputs()
        with patch.dict(os.environ, {extract.PARENT_TRACKING_ENV: str(dwi)}), \
                patch.object(extract, "start_output_tracking") as start:
            status, stdout, stderr = self.run_main(image, roi)
        self.assertEqual(status, 0, stderr)
        self.assertEqual(stderr, "")
        start.assert_not_called()
        self.assertIn("Done:", stdout)
        self.assertIn("Manifest tracking: parent batch", stdout)
        self.assertTrue(Path(extract.getOutfile(roi, image)).is_file())
        log = (dwi / "sub-NR1037_ses-PT3_dwi_analysis.log").read_text()
        self.assertIn("Manifest tracking: parent batch", log)
        self.assertIn("Done:", log)
        self.assertEqual(list(dwi.glob(".*aidamri_manifest.json")), [])

    def test_foreign_or_relative_parent_marker_keeps_standalone_tracking(self):
        dwi, image, roi = self.bids_inputs()
        original_start = extract.start_output_tracking
        for marker in (str(self.folder / "another-dwi"), ".", "1", ""):
            with self.subTest(marker=marker), patch.dict(os.environ, {extract.PARENT_TRACKING_ENV: marker}), \
                    patch.object(extract, "start_output_tracking", wraps=original_start) as start:
                status, stdout, stderr = self.run_main(image, roi)
                self.assertEqual(status, 0, stderr)
                self.assertEqual(stderr, "")
                start.assert_called_once_with(dwi, "dwi", "extraction")
                self.assertNotIn("Manifest tracking: parent batch", stdout)
                self.assertIn("DSI_studio/scan.fz.fa_Anno.csv",
                              self.read_manifest(dwi)["stages"]["extraction"]["created_files"])


if __name__ == "__main__":
    unittest.main()
