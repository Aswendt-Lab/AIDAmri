"""Create a project-level single-image or overlay report from NIfTI suffixes."""

import argparse
import logging
from pathlib import Path

if __package__:
    from . import batch_qc_reports as qc
else:
    import batch_qc_reports as qc

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np


GEOMETRY_ATOL = 1e-4
GEOMETRY_RTOL = 1e-5
MODALITY_CHOICES = ("all", "anat", "dwi", "func", "t2map")


def _nifti_suffix(value):
    if not value or "/" in value or "\\" in value or not value.endswith((".nii", ".nii.gz")):
        raise argparse.ArgumentTypeError("must be a filename suffix ending in .nii or .nii.gz")
    return value


def _opacity(value):
    value = float(value)
    if not np.isfinite(value) or not 0 <= value <= 100:
        raise argparse.ArgumentTypeError("must be a percentage between 0 and 100")
    return value


def _tolerance(value):
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError("must be finite and nonnegative")
    return value


def _load_volume(path):
    image = nib.load(str(path))
    data = image.get_fdata()
    if data.ndim == 4:
        data = data[..., 0]
    if data.ndim != 3:
        raise ValueError(f"Expected 3D or 4D NIfTI data, got shape {image.shape}")
    return image, data


def _check_geometry(image1, image2, atol=GEOMETRY_ATOL, rtol=GEOMETRY_RTOL):
    """Require the same voxel grid; small header rounding differences are allowed."""
    if image1.shape[:3] != image2.shape[:3]:
        raise ValueError(f"Shape mismatch: {image1.shape[:3]} vs {image2.shape[:3]}")
    units1 = image1.header.get_xyzt_units()[0]
    units2 = image2.header.get_xyzt_units()[0]
    if units1 != units2 and "unknown" not in (units1, units2):
        raise ValueError(f"Spatial unit mismatch: {units1} vs {units2}")
    zooms1 = np.asarray(image1.header.get_zooms()[:3])
    zooms2 = np.asarray(image2.header.get_zooms()[:3])
    for zooms in (zooms1, zooms2):
        if not np.all(np.isfinite(zooms)) or np.any(zooms <= 0):
            raise ValueError("Voxel sizes must be finite and positive")
    if not np.allclose(zooms1, zooms2, atol=atol, rtol=rtol):
        raise ValueError(f"Voxel size mismatch: {tuple(zooms1)} vs {tuple(zooms2)}")
    if not all(np.all(np.isfinite(image.affine)) for image in (image1, image2)):
        raise ValueError("Spatial affines must be finite")
    if nib.aff2axcodes(image1.affine) != nib.aff2axcodes(image2.affine):
        raise ValueError("Orientation mismatch between the two NIfTI images")
    if not np.allclose(image1.affine, image2.affine, atol=atol, rtol=rtol):
        raise ValueError(
            f"Spatial affine mismatch (position, orientation or voxel grid; atol={atol}, rtol={rtol})"
        )


def _plot_nifti_image(
    nifti_path, overlay_path, out_dir, project_dir, n_slices, opacity,
    geometry_atol, geometry_rtol, filename_suffix,
):
    image, data = _load_volume(nifti_path)
    overlay_image = overlay_data = None
    if overlay_path is not None:
        overlay_image, overlay_data = _load_volume(overlay_path)
        _check_geometry(image, overlay_image, atol=geometry_atol, rtol=geometry_rtol)

    zooms = image.header.get_zooms()[:3]
    orientations = ["Axial", "Sagittal", "Coronal"]
    slices = qc._slice_indices(data, n_slices)
    vmin, vmax = qc._display_limits(data)
    alpha = opacity / 100
    fig, axes = plt.subplots(3, n_slices, figsize=(3 * n_slices, 9))
    axes = np.asarray(axes).reshape(3, n_slices)
    try:
        for row, orientation in enumerate(orientations):
            for col, index in enumerate(slices[row]):
                ax = axes[row, col]
                img_slice = np.rot90(qc._slice(data, orientation, index))
                extent, x, y = qc._display_geometry(img_slice, orientation, zooms)
                ax.imshow(
                    img_slice, cmap="gray", vmin=vmin, vmax=vmax,
                    extent=extent, origin="lower", aspect="equal",
                )
                if overlay_data is not None and alpha > 0:
                    overlay_slice = np.rot90(qc._slice(overlay_data, orientation, index))
                    foreground = np.isfinite(overlay_slice) & (overlay_slice > 0)
                    overlay = np.ma.masked_where(~foreground, overlay_slice)
                    ax.imshow(
                        overlay, cmap="tab20", alpha=alpha, interpolation="nearest",
                        extent=extent, origin="lower", aspect="equal",
                    )
                    if np.any(foreground) and np.any(~foreground):
                        ax.contour(
                            x, y, foreground, levels=[0.5], colors="yellow",
                            linewidths=0.45, alpha=alpha,
                        )
                ax.set_title(f"{orientation} {index}", fontsize=9)
                ax.axis("off")

        report_title = "Overlay Report" if overlay_path is not None else "Display Report"
        title = f"{report_title}: {Path(nifti_path).name}"
        if overlay_path is not None:
            title += f" + {Path(overlay_path).name}"
        fig.suptitle(title, fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        png_path = Path(out_dir) / qc._safe_png_name(nifti_path, project_dir, filename_suffix)
        fig.savefig(png_path, dpi=120)
    finally:
        plt.close(fig)
    return png_path, image.shape, overlay_image.shape if overlay_image is not None else None, zooms


def build_display_nifti_report(
    project_dir, nifti_in1, nifti_in2=None, opacity=35, n_slices=10,
    geometry_atol=GEOMETRY_ATOL, geometry_rtol=GEOMETRY_RTOL, modality="all",
):
    project_dir = Path(project_dir).expanduser().resolve()
    if not project_dir.is_dir():
        raise ValueError(f"Project directory does not exist or is not a directory: {project_dir}")
    if modality not in MODALITY_CHOICES:
        raise ValueError(f"Unknown modality {modality!r}; choose from {', '.join(MODALITY_CHOICES)}")
    _nifti_suffix(nifti_in1)
    if nifti_in2 is not None:
        _nifti_suffix(nifti_in2)
    opacity = _opacity(opacity)
    geometry_atol = _tolerance(geometry_atol)
    geometry_rtol = _tolerance(geometry_rtol)
    n_slices = qc._positive_int(n_slices)

    report_prefix = "overlay_nifti_report" if nifti_in2 is not None else "display_nifti_report"
    report_stem = f"{report_prefix}_{nifti_in1}"
    if nifti_in2 is not None:
        report_stem += f"_{nifti_in2}"
    out_dir = project_dir / "Report" / ("Overlay" if nifti_in2 is not None else "Display")
    out_dir.mkdir(parents=True, exist_ok=True)
    entries = []
    modality_pattern = "*" if modality == "all" else modality
    for modality_dir in sorted(project_dir.glob(f"sub-*/ses-*/{modality_pattern}")):
        if not modality_dir.is_dir():
            continue
        files1, files2 = [], []
        for path in sorted(modality_dir.rglob("*")):
            if not path.is_file():
                continue
            if path.name.endswith(nifti_in1):
                files1.append(path)
            if nifti_in2 is not None and path.name.endswith(nifti_in2):
                files2.append(path)

        if nifti_in2 is None:
            pairs = [(path, None) for path in files1]
        else:
            if not files1 and not files2:
                continue
            rel_modality = modality_dir.relative_to(project_dir)
            if len(files1) > 1 or len(files2) > 1:
                logging.warning(
                    "Skipping ambiguous combination in %s: %r matches %s; %r matches %s",
                    rel_modality, nifti_in1,
                    [str(path.relative_to(project_dir)) for path in files1],
                    nifti_in2, [str(path.relative_to(project_dir)) for path in files2],
                )
                continue
            if not files1 or not files2:
                logging.warning(
                    "Skipping incomplete combination in %s: %r has %d match(es), %r has %d match(es)",
                    rel_modality, nifti_in1, len(files1), nifti_in2, len(files2),
                )
                continue
            pairs = [(files1[0], files2[0])]

        subject, session, entry_modality = modality_dir.relative_to(project_dir).parts
        for nifti_path, overlay_path in pairs:
            try:
                png_path, shape, overlay_shape, zooms = _plot_nifti_image(
                    nifti_path, overlay_path, out_dir, project_dir, n_slices, opacity,
                    geometry_atol, geometry_rtol, report_stem,
                )
                rel = nifti_path.relative_to(project_dir)
                info = [
                    ("NIfTI 1", rel), ("Modality", entry_modality),
                    ("Dimensions", shape),
                    ("Voxel size", tuple(round(float(z), 4) for z in zooms)),
                ]
                image_alt = nifti_path.name
                if overlay_path is not None:
                    info.extend([
                        ("NIfTI 2", overlay_path.relative_to(project_dir)),
                        ("Overlay dimensions", overlay_shape),
                        ("Overlay opacity", f"{opacity:g}%"),
                    ])
                    image_alt += f" + {overlay_path.name}"
                entries.append({
                    "subject": subject, "session": session, "modality": entry_modality,
                    "report_img_path": png_path.name, "image_alt": image_alt, "info": info,
                })
            except Exception as exc:
                logging.warning(
                    "Skipping Display Report for %s%s: %s", nifti_path,
                    f" + {overlay_path}" if overlay_path is not None else "", exc,
                )

    if not entries:
        return None, 0
    report_title = "Overlay Report" if nifti_in2 is not None else "Display Report"
    title = f"{report_title}: {nifti_in1}"
    if nifti_in2 is not None:
        title += f" + {nifti_in2}"
    html_path = qc._write_report(
        entries, out_dir, title, f"{report_stem}.html",
    )
    return html_path, len(entries)


def _build_argument_parser():
    parser = argparse.ArgumentParser(
        description="Create a project-level NIfTI display or overlay report using filename suffixes."
    )
    parser.add_argument(
        "-i", "--input", dest="project_dir", required=True, type=Path,
        help="Project directory containing sub-*/ses-*/<modality> folders.",
    )
    parser.add_argument(
        "-1", "--nifti_in1", required=True, type=_nifti_suffix,
        help="Base-image filename suffix, e.g. fz.fa.nii.gz; searched recursively per modality.",
    )
    parser.add_argument(
        "-2", "--nifti_in2", type=_nifti_suffix,
        help="Optional overlay suffix; requires exactly one match per input in the same modality.",
    )
    parser.add_argument(
        "-m", "--modality", choices=MODALITY_CHOICES, default="all",
        help="Modality folders to search, including their subfolders (default: all).",
    )
    parser.add_argument(
        "-o", "--opacity", type=_opacity, default=None,
        help="Opacity of NIfTI 2 in percent, 0–100 (default: 35); requires --nifti_in2.",
    )
    parser.add_argument(
        "-n", "--n-slices", type=qc._positive_int, default=10,
        help="Number of slices per orientation (default: 10).",
    )
    parser.add_argument(
        "--geometry-atol", type=_tolerance, default=GEOMETRY_ATOL,
        help="Absolute tolerance for voxel sizes and spatial affines (default: 1e-4, header units).",
    )
    parser.add_argument(
        "--geometry-rtol", type=_tolerance, default=GEOMETRY_RTOL,
        help="Relative tolerance for voxel sizes and spatial affines (default: 1e-5).",
    )
    return parser


def main(argv=None):
    parser = _build_argument_parser()
    args = parser.parse_args(argv)
    project_dir = args.project_dir.expanduser()
    if not project_dir.is_dir():
        parser.error(f"project directory does not exist or is not a directory: {project_dir}")
    if args.opacity is not None and args.nifti_in2 is None:
        parser.error("--opacity requires --nifti_in2")
    opacity = 35 if args.opacity is None else args.opacity
    html_path, count = build_display_nifti_report(
        project_dir, args.nifti_in1, nifti_in2=args.nifti_in2, opacity=opacity,
        n_slices=args.n_slices, geometry_atol=args.geometry_atol,
        geometry_rtol=args.geometry_rtol, modality=args.modality,
    )
    if html_path:
        print(f"Display report written to {html_path} ({count} image(s))")
    else:
        print("Display report skipped: no matching files or compatible pairs found.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
