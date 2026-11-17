"""Check converted input selection and the shared report presentation."""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from helper_tools import batch_qc_reports as reports
from matplotlib.axes import Axes
from matplotlib.figure import Figure


class ConversionReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.project = Path(self.temp.name)

    def save(self, relative, data=None):
        path = self.project / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if data is None:
            data = np.ones((3, 4, 5))
        nib.save(nib.Nifti1Image(data, np.diag([0.2, 0.3, 0.4, 1])), path)
        return path

    def test_selection_metadata_shared_viewer_and_unique_images(self):
        name = "sub-01_ses-01_scan.nii.gz"
        for folder in ("anat", "func", "dwi", "dwi/derived"):
            self.save(f"sub-01/ses-01/{folder}/{name}")
        self.save("sub-02/ses-02/func/sub-02_ses-02_bold.nii", np.ones((3, 4, 5, 2)))
        for folder, filename in (("t2map", name), ("anat/deep", name),
                                 ("dwi/deep/deeper", name), ("anat", "other.nii")):
            self.save(f"sub-01/ses-01/{folder}/{filename}")
        html_path, count = reports.build_conversion_qc_report(self.project, n_slices=1)
        self.assertEqual(count, 5)
        self.assertEqual(html_path, self.project / "Report/Convert2Nifti/Convert2Nifti_Report.html")
        html = html_path.read_text()
        for text in ("image-viewer", "viewer-zoom-in", "subjectDropdown", "sessionDropdown",
                     "modalityDropdown", "# Volumes:</b> 2", "(0.2, 0.3, 0.4)"):
            self.assertIn(text, html)
        self.assertNotIn("max-width: 1200px", html)
        pngs = list(html_path.parent.glob("*.png"))
        self.assertEqual(len(pngs), count)
        for png in pngs:
            self.assertIn(png.name, html)
            self.assertTrue(png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

    def test_default_seven_slices_first_volume_and_geometry(self):
        first = np.arange(60, dtype=float).reshape(3, 4, 5)
        self.save("sub-01/ses-01/func/sub-01_ses-01_bold.nii.gz",
                  np.stack([first, first + 1000], axis=-1))
        images = []
        imshow = Axes.imshow

        def capture(ax, data, *args, **kwargs):
            images.append((np.asarray(data).copy(), kwargs))
            return imshow(ax, data, *args, **kwargs)

        # Keep this test fast; the other test verifies actual PNG encoding.
        with patch.object(Axes, "imshow", capture), patch.object(Figure, "savefig"):
            _, count = reports.build_conversion_qc_report(self.project)
        self.assertEqual(count, 1)
        self.assertEqual(len(images), 21)
        np.testing.assert_array_equal(images[0][0], np.rot90(first[:, :, 0]))
        np.testing.assert_allclose(images[0][1]["extent"], (0, 0.6, 0, 1.2))
        self.assertTrue(all(image.max() < 1000 for image, _ in images))

    def test_lip_labels_match_existing_rotated_views(self):
        data = np.arange(60, dtype=float).reshape(3, 4, 5)
        affine = np.array([[-0.2, 0, 0, 0], [0, 0, -0.4, 0],
                           [0, -0.3, 0, 0], [0, 0, 0, 1]])
        path = self.project / "lip.nii.gz"
        nib.save(nib.Nifti1Image(data, affine), path)
        before = path.read_bytes()
        with patch.object(Figure, "savefig"), \
                patch.object(Axes, "imshow") as show, \
                patch.object(Axes, "annotate") as annotate:
            reports._plot_conversion_image(path, self.project, self.project, 1)
        expected_planes = (data[:, :, 0], data[0, :, :], data[:, 0, :])
        for call, plane in zip(show.call_args_list, expected_planes):
            np.testing.assert_array_equal(call.args[0], np.rot90(plane))
        actual_labels = [(call.args[0], call.kwargs["xy"]) for call in annotate.call_args_list]
        expected_labels = [
            ("R", (0, 0.5)), ("L", (1, 0.5)), ("S", (0.5, 1)), ("I", (0.5, 0)),
            ("S", (0, 0.5)), ("I", (1, 0.5)), ("A", (0.5, 1)), ("P", (0.5, 0)),
            ("R", (0, 0.5)), ("L", (1, 0.5)), ("A", (0.5, 1)), ("P", (0.5, 0)),
        ]
        self.assertEqual(actual_labels, expected_labels)
        self.assertEqual(path.read_bytes(), before)

    def test_empty_and_corrupt_inputs(self):
        self.assertEqual(reports.build_conversion_qc_report(self.project), (None, 0))
        path = self.project / "sub-01/ses-01/anat/sub-broken.nii.gz"
        path.parent.mkdir(parents=True)
        path.write_bytes(b"not a nifti")
        with self.assertLogs(level="WARNING"):
            self.assertEqual(reports.build_conversion_qc_report(self.project), (None, 0))
        self.save("sub-01/ses-01/anat/sub-valid.nii.gz")
        with self.assertLogs(level="WARNING"):
            _, count = reports.build_conversion_qc_report(self.project, n_slices=1)
        self.assertEqual(count, 1)


if __name__ == "__main__":
    unittest.main()
