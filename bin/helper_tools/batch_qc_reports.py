import argparse
import html
import logging
import os
from calendar import month_name
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

REPORT_TIMEZONE = ZoneInfo("Europe/Berlin")


def _as_3d(data):
    data = np.asarray(data)
    if data.ndim == 4:
        data = data[..., 0]
    data = np.squeeze(data)
    if data.ndim != 3:
        raise ValueError(f"Expected 3D or 4D NIfTI data, got shape {data.shape}")
    return data


def _load_3d(nifti_path):
    img = nib.load(str(nifti_path))
    zooms = img.header.get_zooms()[:3]
    return _as_3d(img.get_fdata()), img.shape, zooms


def _slice_indices(data, n_slices):
    return [
        np.linspace(0, data.shape[2] - 1, n_slices, dtype=int),
        np.linspace(0, data.shape[0] - 1, n_slices, dtype=int),
        np.linspace(0, data.shape[1] - 1, n_slices, dtype=int),
    ]


def _slice(data, orientation, index):
    if orientation == "Axial":
        return data[:, :, index]
    if orientation == "Sagittal":
        return data[index, :, :]
    return data[:, index, :]


def _display_geometry(img_slice, orientation, zooms):
    if orientation == "Axial":
        row_spacing, col_spacing = zooms[1], zooms[0]
    elif orientation == "Sagittal":
        row_spacing, col_spacing = zooms[2], zooms[1]
    else:
        row_spacing, col_spacing = zooms[2], zooms[0]

    rows, cols = img_slice.shape
    extent = (0, cols * col_spacing, 0, rows * row_spacing)
    x = np.linspace(0, cols * col_spacing, cols)
    y = np.linspace(0, rows * row_spacing, rows)
    return extent, x, y


def _display_limits(data):
    finite = data[np.isfinite(data)]
    nonzero = finite[finite != 0]
    values = nonzero if nonzero.size else finite
    if values.size == 0:
        return 0, 1
    low, high = np.percentile(values, [1, 99])
    if low == high:
        high = low + 1
    return low, high


def _safe_png_name(nifti_path, project_dir, suffix):
    rel = Path(nifti_path).resolve().relative_to(Path(project_dir).resolve())
    name = str(rel).replace(os.sep, "__")
    name = name.replace(".nii.gz", "").replace(".nii", "")
    return f"{name}_{suffix}.png"


def _entry_metadata(nifti_path, project_dir):
    rel = Path(nifti_path).resolve().relative_to(Path(project_dir).resolve())
    parts = rel.parts
    subject = next((part for part in parts if part.startswith("sub-")), "sub-unknown")
    session = next((part for part in parts if part.startswith("ses-")), "ses-unknown")
    modality = "unknown"
    for candidate in ("anat", "dwi", "func", "t2map"):
        if candidate in parts:
            modality = candidate
            break
    return rel, subject, session, modality


def _plot_bet_image(bet_path, out_dir, project_dir, n_slices):
    data, shape, zooms = _load_3d(bet_path)
    mask_path = Path(str(bet_path).replace(".nii.gz", "_mask.nii.gz"))
    mask = None
    if mask_path.exists():
        try:
            mask, _, _ = _load_3d(mask_path)
            if mask.shape != data.shape:
                logging.warning("Skipping BET mask overlay with mismatched shape: %s", mask_path)
                mask = None
        except Exception as exc:
            logging.warning("Could not load BET mask %s: %s", mask_path, exc)

    orientations = ["Axial", "Sagittal", "Coronal"]
    slices = _slice_indices(data, n_slices)
    vmin, vmax = _display_limits(data)
    fig, axes = plt.subplots(3, n_slices, figsize=(3 * n_slices, 9))
    axes = np.asarray(axes).reshape(3, n_slices)

    for row, orientation in enumerate(orientations):
        for col, index in enumerate(slices[row]):
            ax = axes[row, col]
            img_slice = np.rot90(_slice(data, orientation, index))
            extent, x, y = _display_geometry(img_slice, orientation, zooms)
            ax.imshow(
                img_slice,
                cmap="gray",
                vmin=vmin,
                vmax=vmax,
                extent=extent,
                origin="lower",
                aspect="equal",
            )
            if mask is not None:
                mask_slice = np.rot90(_slice(mask, orientation, index) > 0)
                if np.any(mask_slice) and np.any(~mask_slice):
                    ax.contour(x, y, mask_slice, levels=[0.5], colors="red", linewidths=0.7)
            ax.set_title(f"{orientation} {index}", fontsize=9)
            ax.axis("off")

    fig.suptitle(f"BET Report: {Path(bet_path).name}", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    png_path = Path(out_dir) / _safe_png_name(bet_path, project_dir, "bet_report")
    fig.savefig(png_path, dpi=120)
    plt.close(fig)
    return png_path, shape, zooms, mask_path if mask_path.exists() else None


def _plot_registration_overlay(bet_path, anno_path, out_dir, project_dir, n_slices):
    bet_data, shape, zooms = _load_3d(bet_path)
    anno_data, anno_shape, _ = _load_3d(anno_path)
    if bet_data.shape != anno_data.shape:
        raise ValueError(
            f"Shape mismatch for overlay: BET {bet_data.shape}, annotation {anno_data.shape}"
        )

    orientations = ["Axial", "Sagittal", "Coronal"]
    slices = _slice_indices(bet_data, n_slices)
    vmin, vmax = _display_limits(bet_data)
    fig, axes = plt.subplots(3, n_slices, figsize=(3 * n_slices, 9))
    axes = np.asarray(axes).reshape(3, n_slices)

    for row, orientation in enumerate(orientations):
        for col, index in enumerate(slices[row]):
            ax = axes[row, col]
            bet_slice = np.rot90(_slice(bet_data, orientation, index))
            anno_slice = np.rot90(_slice(anno_data, orientation, index))
            anno_overlay = np.ma.masked_where(anno_slice <= 0, anno_slice)
            extent, x, y = _display_geometry(bet_slice, orientation, zooms)
            ax.imshow(
                bet_slice,
                cmap="gray",
                vmin=vmin,
                vmax=vmax,
                extent=extent,
                origin="lower",
                aspect="equal",
            )
            ax.imshow(
                anno_overlay,
                cmap="tab20",
                alpha=0.35,
                interpolation="nearest",
                extent=extent,
                origin="lower",
                aspect="equal",
            )
            if np.any(anno_slice > 0) and np.any(anno_slice <= 0):
                ax.contour(x, y, anno_slice > 0, levels=[0.5], colors="yellow", linewidths=0.45)
            ax.set_title(f"{orientation} {index}", fontsize=9)
            ax.axis("off")

    fig.suptitle(f"Registration Report: {Path(bet_path).name} + {Path(anno_path).name}", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    png_path = Path(out_dir) / _safe_png_name(anno_path, project_dir, "registration_report")
    fig.savefig(png_path, dpi=120)
    plt.close(fig)
    return png_path, shape, anno_shape, zooms


def _format_parameter_value(value):
    if isinstance(value, (list, tuple)):
        return " ".join(str(item) for item in value)
    return str(value)


def _report_timestamp():
    timestamp = datetime.now(REPORT_TIMEZONE)
    return (
        f"{timestamp.day:02d} {month_name[timestamp.month]} {timestamp.year} "
        f"{timestamp:%H:%M:%S} {timestamp.tzname()}"
    )


def _write_report(entries, out_dir, title, report_name, custom_parameters=None):
    html_path = Path(out_dir) / report_name
    subjects = sorted({entry["subject"] for entry in entries})
    sessions = sorted({entry["session"] for entry in entries})
    modalities = sorted({entry["modality"] for entry in entries})
    generated_at = _report_timestamp()

    with open(html_path, "w") as f:
        f.write(f"<html><head><title>{html.escape(title)}</title>\n")
        f.write(
            """
            <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .custom-parameters { background: #f3f6f8; border: 1px solid #ccd5db; border-radius: 4px; margin-bottom: 24px; padding: 16px 20px; }
            .custom-parameters h2 { margin: 0 0 12px; }
            .custom-parameters table { border-collapse: collapse; }
            .custom-parameters th { padding: 3px 24px 3px 0; text-align: left; vertical-align: top; }
            .custom-parameters td { padding: 3px 0; }
            .report-generated { color: #555; font-size: 0.95em; margin: -10px 0 24px; }
            .report-entry { margin-bottom: 40px; }
            .report-info { font-size: 1.0em; margin-bottom: 8px; line-height: 1.5; }
            .report-img { width: 100%; max-width: 1200px; border: 1px solid #ccc; }
            #report-dropdown-bar { position: fixed; top: 0; left: 0; width: 100%; background: #f9f9f9; border-bottom: 1px solid #ccc; z-index: 1000; padding: 12px 40px; box-sizing: border-box; }
            </style>
            <script>
            function filterreport() {
                var subj = document.getElementById('subjectDropdown').value;
                var sess = document.getElementById('sessionDropdown').value;
                var mod = document.getElementById('modalityDropdown').value;
                var entries = document.getElementsByClassName('report-entry');
                for (var i = 0; i < entries.length; i++) {
                    var entry = entries[i];
                    var show = true;
                    if (subj !== 'all' && entry.getAttribute('data-subject') !== subj) show = false;
                    if (sess !== 'all' && entry.getAttribute('data-session') !== sess) show = false;
                    if (mod !== 'all' && entry.getAttribute('data-modality') !== mod) show = false;
                    entry.style.display = show ? '' : 'none';
                }
            }
            </script>
            """
        )
        f.write("</head><body>\n")
        f.write("<div id='report-dropdown-bar'>\n")
        for label, element_id, values in (
            ("Subject", "subjectDropdown", subjects),
            ("Session", "sessionDropdown", sessions),
            ("Modality", "modalityDropdown", modalities),
        ):
            f.write(f"<label style='margin-right:20px;'>{label}: ")
            f.write(f"<select id='{element_id}' onchange='filterreport()'>")
            f.write("<option value='all'>All</option>")
            for value in values:
                escaped_value = html.escape(value)
                f.write(f"<option value='{escaped_value}'>{escaped_value}</option>")
            f.write("</select></label>\n")
        f.write("</div><div style='height:60px;'></div>\n")
        f.write(f"<h1>{html.escape(title)}</h1>\n")
        f.write(f"<p class='report-generated'><b>Created:</b> {html.escape(generated_at)}</p>\n")
        if custom_parameters:
            f.write("<section class='custom-parameters'>\n")
            f.write("<h2>Custom parameters</h2>\n<table>\n")
            for name, value in custom_parameters:
                f.write(
                    "<tr>"
                    f"<th>{html.escape(str(name))}</th>"
                    f"<td>{html.escape(_format_parameter_value(value))}</td>"
                    "</tr>\n"
                )
            f.write("</table>\n</section>\n")

        for entry in entries:
            f.write(
                "<div class='report-entry' "
                f"data-subject='{html.escape(entry['subject'])}' "
                f"data-session='{html.escape(entry['session'])}' "
                f"data-modality='{html.escape(entry['modality'])}'>\n"
            )
            f.write("<div class='report-info'>")
            for label, value in entry["info"]:
                f.write(f"<b>{html.escape(label)}:</b> {html.escape(str(value))} &nbsp; ")
            f.write("</div>\n")
            f.write(
                f"<img class='report-img' src='{html.escape(entry['report_img_path'])}' "
                f"alt='{html.escape(entry['image_alt'])}'>\n"
            )
            f.write("</div>\n")
        f.write("</body></html>\n")
    return html_path


def build_bet_qc_report(project_dir, n_slices=10, custom_parameters=None):
    project_dir = Path(project_dir)
    out_dir = project_dir / "Report" / "BET"
    out_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    bet_files = sorted(project_dir.glob("sub-*/ses-*/*/*Bet.nii.gz"))
    for bet_path in bet_files:
        if bet_path.name.endswith("_mask.nii.gz"):
            continue
        try:
            png_path, shape, zooms, mask_path = _plot_bet_image(bet_path, out_dir, project_dir, n_slices)
            rel, subject, session, modality = _entry_metadata(bet_path, project_dir)
            info = [
                ("File", rel),
                ("Modality", modality),
                ("Dimensions", shape),
                ("Voxel size", tuple(round(float(z), 4) for z in zooms)),
                ("Mask", mask_path.name if mask_path else "not found"),
            ]
            entries.append(
                {
                    "subject": subject,
                    "session": session,
                    "modality": modality,
                    "report_img_path": png_path.name,
                    "image_alt": bet_path.name,
                    "info": info,
                }
            )
        except Exception as exc:
            logging.warning("Could not create BET Report for %s: %s", bet_path, exc)

    if not entries:
        return None, 0
    html_path = _write_report(
        entries,
        out_dir,
        "BET Report",
        "bet_report.html",
        custom_parameters=custom_parameters,
    )
    return html_path, len(entries)


def build_registration_qc_report(project_dir, n_slices=10, custom_parameters=None):
    project_dir = Path(project_dir)
    out_dir = project_dir / "Report" / "Registration"
    out_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    suffix = "_AnnoSplit_parental.nii.gz"
    anno_files = sorted(project_dir.glob(f"sub-*/ses-*/*/*{suffix}"))
    for anno_path in anno_files:
        bet_name = anno_path.name[: -len(suffix)] + ".nii.gz"
        bet_path = anno_path.with_name(bet_name)
        if not bet_path.exists():
            logging.warning("Skipping registration Report without matching BET file: %s", anno_path)
            continue
        try:
            png_path, bet_shape, anno_shape, zooms = _plot_registration_overlay(
                bet_path, anno_path, out_dir, project_dir, n_slices
            )
            rel_bet, subject, session, modality = _entry_metadata(bet_path, project_dir)
            rel_anno = anno_path.resolve().relative_to(project_dir.resolve())
            info = [
                ("BET", rel_bet),
                ("Annotation", rel_anno),
                ("Modality", modality),
                ("BET dimensions", bet_shape),
                ("Annotation dimensions", anno_shape),
                ("Voxel size", tuple(round(float(z), 4) for z in zooms)),
            ]
            entries.append(
                {
                    "subject": subject,
                    "session": session,
                    "modality": modality,
                    "report_img_path": png_path.name,
                    "image_alt": f"{bet_path.name} + {anno_path.name}",
                    "info": info,
                }
            )
        except Exception as exc:
            logging.warning("Could not create registration Report for %s: %s", anno_path, exc)

    if not entries:
        return None, 0
    html_path = _write_report(
        entries,
        out_dir,
        "Registration Report: BET + AnnoSplit_parental",
        "registration_report.html",
        custom_parameters=custom_parameters,
    )
    return html_path, len(entries)


def _positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return value


def _custom_parameter(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError("must use NAME=VALUE format")

    name, parameter_value = value.split("=", 1)
    name = name.strip()
    if not name:
        raise argparse.ArgumentTypeError("parameter name must not be empty")
    if not name.startswith("-"):
        name = f"--{name}"
    return name, parameter_value


def _build_argument_parser():
    parser = argparse.ArgumentParser(
        description="Create project-level BET and registration QC reports."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="project_dir",
        required=True,
        type=Path,
        help="Processed project directory containing sub-*/ses-* folders.",
    )
    parser.add_argument(
        "--report",
        choices=("all", "bet", "registration"),
        default="all",
        help="Report to create (default: all).",
    )
    parser.add_argument(
        "--n-slices",
        type=_positive_int,
        default=10,
        help="Number of slices per orientation (default: 10).",
    )
    parser.add_argument(
        "--custom-parameter",
        dest="custom_parameters",
        action="append",
        type=_custom_parameter,
        default=[],
        metavar="NAME=VALUE",
        help=(
            "Parameter to list in the HTML report, for example "
            "--custom-parameter t2-frac=0.1. May be repeated."
        ),
    )
    return parser


def main(argv=None):
    parser = _build_argument_parser()
    args = parser.parse_args(argv)
    project_dir = args.project_dir.expanduser()

    if not project_dir.is_dir():
        parser.error(f"project directory does not exist or is not a directory: {project_dir}")

    report_builders = []
    if args.report in ("all", "bet"):
        report_builders.append(("BET", build_bet_qc_report))
    if args.report in ("all", "registration"):
        report_builders.append(("Registration", build_registration_qc_report))

    for report_label, report_builder in report_builders:
        html_path, count = report_builder(
            project_dir,
            n_slices=args.n_slices,
            custom_parameters=args.custom_parameters or None,
        )
        if html_path:
            print(f"{report_label} report written to {html_path} ({count} image(s))")
        else:
            print(f"{report_label} report skipped: no matching files found.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
