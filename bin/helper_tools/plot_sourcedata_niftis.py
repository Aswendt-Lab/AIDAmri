import os
import argparse
from calendar import month_name
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

REPORT_TIMEZONE = ZoneInfo("Europe/Berlin")

def _display_geometry(img_slice, orientation, zooms):
    if orientation == "Axial":
        row_spacing, col_spacing = zooms[1], zooms[0]
    elif orientation == "Sagittal":
        row_spacing, col_spacing = zooms[2], zooms[1]
    else:
        row_spacing, col_spacing = zooms[2], zooms[0]

    rows, cols = img_slice.shape
    extent = (0, cols * col_spacing, 0, rows * row_spacing)
    return extent

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

def _voxel_size(img):
    return tuple(round(float(z), 4) for z in img.header.get_zooms()[:3])

def plot_nifti_slices(nifti_path, out_dir, n_slices=10):
    img = nib.load(nifti_path)
    zooms = img.header.get_zooms()[:3]
    # Attempt to load data, if end of file error occurs, skip this file and create a placeholder image
    try:
        data = img.get_fdata()
    except EOFError as e:
        print(f"Error loading {nifti_path}: {e}")
        # Create a placeholder image
        data = np.zeros((10, 10, 10))  # Placeholder for empty data
        plt.imshow(data[:, :, 0], cmap='gray')  # Show a single slice
        plt.title(f"Placeholder for {os.path.basename(nifti_path)}")
        plt.axis('off')
        out_path = os.path.join(out_dir, os.path.basename(nifti_path).replace('.nii', '').replace('.gz', '') + '_report.png')
        plt.savefig(out_path)
        plt.close()
        return out_path
    fname = os.path.basename(nifti_path)
    fig, axes = plt.subplots(3, n_slices, figsize=(3*n_slices, 9))
    orientations = ['Axial', 'Sagittal', 'Coronal']
    slices = [
        np.linspace(0, data.shape[2]-1, n_slices, dtype=int),  # Axial (z)
        np.linspace(0, data.shape[0]-1, n_slices, dtype=int),  # Sagittal (x)
        np.linspace(0, data.shape[1]-1, n_slices, dtype=int)   # Coronal (y)
    ]
    # Get orientation labels based on NIfTI header
    affine = img.affine
    # Use nibabel's orientation info to get axis codes (e.g., 'L', 'R', 'A', 'P', 'S', 'I')
    axis_codes = nib.orientations.aff2axcodes(affine)
    # Map orientations to axes for each view
    # Axial: slices along z, show x/y axes
    # Sagittal: slices along x, show y/z axes
    # Coronal: slices along y, show x/z axes
    view_axes = [
        (0, 1, 2),  # Axial: x/y, slice z
        (1, 2, 0),  # Sagittal: y/z, slice x
        (0, 2, 1)   # Coronal: x/z, slice y
    ]
    # If data is 4D, use the first volume
    if data.ndim == 4:
        data = data[..., 0]
    vmin, vmax = _display_limits(data)
    for i, (ori, slcs) in enumerate(zip(orientations, slices)):
        ax1, ax2, _ = view_axes[i]
        # axis_codes gives the positive direction for each axis; get negative direction by swapping L<->R, A<->P, S<->I
        def opposite(code):
            opposites = {'L': 'R', 'R': 'L', 'A': 'P', 'P': 'A', 'S': 'I', 'I': 'S'}
            return opposites.get(code, code)
        left, right = axis_codes[ax1], opposite(axis_codes[ax1])
        top, bottom = axis_codes[ax2], opposite(axis_codes[ax2])
        for j, sl in enumerate(slcs):
            # Select slice based on orientation
            if ori == 'Axial':
                img_slice = data[:, :, sl]
            elif ori == 'Sagittal':
                img_slice = data[sl, :, :]
            else:  # Coronal
                img_slice = data[:, sl, :]
            # Ensure img_slice is 2D
            img_slice = np.rot90(np.squeeze(img_slice))
            extent = _display_geometry(img_slice, ori, zooms)
            axes[i, j].imshow(
                img_slice,
                cmap='gray',
                vmin=vmin,
                vmax=vmax,
                extent=extent,
                origin='lower',
                aspect='equal'
            )
            axes[i, j].set_title(f"{ori} slice {sl}", fontsize=12)
            axes[i, j].axis('off')
            # Add orientation labels from axis_codes
            axes[i, j].annotate(left, xy=(0, 0.5), xycoords='axes fraction',
                                va='center', ha='left', fontsize=12, color='lime')
            axes[i, j].annotate(right, xy=(1, 0.5), xycoords='axes fraction',
                                va='center', ha='right', fontsize=12, color='lime')
            axes[i, j].annotate(top, xy=(0.5, 1), xycoords='axes fraction',
                                va='top', ha='center', fontsize=12, color='lime')
            axes[i, j].annotate(bottom, xy=(0.5, 0), xycoords='axes fraction',
                                va='bottom', ha='center', fontsize=12, color='lime')
    plt.suptitle(f"Slices: {fname}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, fname.replace('.nii', '').replace('.gz', '') + '_report.png')
    plt.savefig(out_path)
    plt.close(fig)
    return out_path

def process_subject(subject_dir, out_dir, n_slices=10):
    report_entries = []
    for session in os.listdir(subject_dir):
        session_dir = os.path.join(subject_dir, session)
        if not os.path.isdir(session_dir):
            continue
        for modality in ['dwi', 'anat', 'func']:
            modality_dir = os.path.join(session_dir, modality)
            if not os.path.isdir(modality_dir):
                continue
            for fname in os.listdir(modality_dir):
                if (fname.endswith('.nii') or fname.endswith('.nii.gz')) and fname.startswith('sub-'):
                    nifti_path = os.path.join(modality_dir, fname)
                    img = nib.load(nifti_path)
                    shape = img.shape
                    n_vols = shape[3] if (modality in ['dwi', 'func'] and len(shape) > 3) else 1
                    voxel_size = _voxel_size(img)
                    report_img_path = plot_nifti_slices(nifti_path, out_dir, n_slices)
                    report_entries.append({
                        'filename': fname,
                        'modality': modality,
                        'dimensions': shape,
                        'voxel_size': voxel_size,
                        'n_volumes': n_vols,
                        'report_img_path': os.path.basename(report_img_path)
                    })
            if modality == 'dwi':
                # Plot niftis in subfolders
                for subfolder in os.listdir(modality_dir):
                    subfolder_path = os.path.join(modality_dir, subfolder)
                    if os.path.isdir(subfolder_path):
                        for fname in os.listdir(subfolder_path):
                            if (fname.endswith('.nii') or fname.endswith('.nii.gz')) and fname.startswith('sub-'):
                                try:
                                    nifti_path = os.path.join(subfolder_path, fname)
                                    img = nib.load(nifti_path)
                                    shape = img.shape
                                    n_vols = shape[3] if (modality in ['dwi'] and len(shape) > 3) else 1
                                    voxel_size = _voxel_size(img)
                                    report_img_path = plot_nifti_slices(nifti_path, out_dir, n_slices)
                                    report_entries.append({
                                    'filename': fname,
                                    'modality': modality,
                                    'dimensions': shape,
                                    'voxel_size': voxel_size,
                                    'n_volumes': n_vols,
                                    'report_img_path': os.path.basename(report_img_path)
                                })
                                except Exception as e:
                                    print(f"Error processing {nifti_path}: {e}")
    # Sort entries by modality and filename
    report_entries.sort(key=lambda x: (x['modality'], x['filename']))
    return report_entries

def _report_timestamp():
    timestamp = datetime.now(REPORT_TIMEZONE)
    return (
        f"{timestamp.day:02d} {month_name[timestamp.month]} {timestamp.year} "
        f"{timestamp:%H:%M:%S} {timestamp.tzname()}"
    )

def write_html_report(report_entries, out_dir):
    # Try to extract subject and session IDs from the first entry's filename path
    subject_id = session_id = "unknown"
    if report_entries:
        first_entry = report_entries[0]['filename']
        parts = first_entry.split('_')
        if len(parts) >= 2:
            subject_id = parts[0].replace('sub-', '')
            session_id = parts[1].replace('ses-', '')
    # Compose report file name and title
    report_fname = "Convert2Nifti_Report.html"
    html_path = os.path.join(out_dir, report_fname)
    report_title = "Convert2Nifti_Report"
    generated_at = _report_timestamp()
    subjects = sorted(set(
        entry['filename'].split('_')[0].replace('sub-', '') for entry in report_entries
    ))
    sessions = sorted(set(
        entry['filename'].split('_')[1].replace('ses-', '') for entry in report_entries
    ))
    modalities = sorted(set(entry['modality'] for entry in report_entries))

    with open(html_path, "w") as f:
        f.write(f"<html><head><title>{report_title}</title>\n")
        f.write("""
        <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .report-entry { margin-bottom: 40px; }
        .report-info { font-size: 1.1em; margin-bottom: 8px; }
        .report-img { width: 100%; max-width: 1200px; border: 1px solid #ccc; }
        .report-generated { color: #555; font-size: 0.95em; margin: -10px 0 24px; }
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
        """)
        f.write("</head><body>\n")
        f.write("""
        <div id='report-dropdown-bar' style='position:fixed;top:0;left:0;width:100%;background:#f9f9f9;border-bottom:1px solid #ccc;z-index:1000;padding:12px 0;'>
            <label style='margin-right:20px;'>Subject:
            <select id='subjectDropdown' onchange='filterreport()'>
                <option value='all'>All</option>
        """)
        for s in subjects:
            f.write(f"<option value='{s}'>{s}</option>\n")
        f.write("""
            </select>
            </label>
            <label style='margin-right:20px;'>Session:
            <select id='sessionDropdown' onchange='filterreport()'>
                <option value='all'>All</option>
        """)
        for s in sessions:
            f.write(f"<option value='{s}'>{s}</option>\n")
        f.write("""
            </select>
            </label>
            <label style='margin-right:20px;'>Modality:
            <select id='modalityDropdown' onchange='filterreport()'>
                <option value='all'>All</option>
        """)
        for m in modalities:
            f.write(f"<option value='{m}'>{m}</option>\n")
        f.write("""
            </select>
            </label>
        </div>
        <div style='height:60px;'></div>
        """)
        f.write(f"<h1>{report_title}</h1>\n")
        f.write(f"<p class='report-generated'><b>Created:</b> {generated_at}</p>\n")
        for entry in report_entries:
            # Extract subject and session from filename
            parts = entry['filename'].split('_')
            subj = parts[0].replace('sub-', '') if len(parts) > 0 else "unknown"
            sess = parts[1].replace('ses-', '') if len(parts) > 1 else "unknown"
            f.write(f"<div class='report-entry' data-subject='{subj}' data-session='{sess}' data-modality='{entry['modality']}'>\n")
            f.write(
                f"<div class='report-info'><b>File Name:</b> {entry['filename']} &nbsp; "
                f"<b>Modality:</b> {entry['modality']} &nbsp; "
                f"<b>Dimensions:</b> {entry['dimensions']} &nbsp; "
                f"<b>Voxel size:</b> {entry['voxel_size']} &nbsp; "
                f"<b># Volumes:</b> {entry['n_volumes']}</div>\n"
            )
            f.write(f"<img class='report-img' src='{entry['report_img_path']}' alt='{entry['filename']}'>\n")
            f.write("</div>\n")
        f.write("</body></html>\n")
    print(f"HTML report written to {html_path}")

def main():
    parser = argparse.ArgumentParser(description="Create plots and HTML report for NIfTI files of a subject.")
    parser.add_argument("subject_dir", help="Path to sub-<subject_id>/ses-<session_id>/ directory")
    parser.add_argument("out_dir", help="Directory to save plots and HTML report")
    parser.add_argument("--n_slices", type=int, default=10, help="Number of slices per orientation")
    args = parser.parse_args()
    report_entries = process_subject(args.subject_dir, args.out_dir, args.n_slices)
    write_html_report(report_entries, args.out_dir)

if __name__ == "__main__":
    main()
