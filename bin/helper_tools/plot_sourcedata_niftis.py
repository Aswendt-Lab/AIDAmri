import os
import argparse
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

def plot_nifti_slices(nifti_path, out_dir, n_slices=10):
    img = nib.load(nifti_path)
    data = img.get_fdata()
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
    for i, (ori, slcs) in enumerate(zip(orientations, slices)):
        ax1, ax2, _ = view_axes[i]
        left, right = axis_codes[ax1], axis_codes[ax1][::-1]
        top, bottom = axis_codes[ax2], axis_codes[ax2][::-1]
        for j, sl in enumerate(slcs):
            # Select slice based on orientation
            if ori == 'Axial':
                img_slice = data[:, :, sl]
            elif ori == 'Sagittal':
                img_slice = data[sl, :, :]
            else:  # Coronal
                img_slice = data[:, sl, :]
            # Ensure img_slice is 2D
            img_slice = np.squeeze(img_slice)
            axes[i, j].imshow(np.rot90(img_slice), cmap='gray')
            axes[i, j].set_title(f"{ori} slice {sl}", fontsize=12)
            axes[i, j].axis('off')
            # Add orientation labels from axis_codes
            axes[i, j].annotate(left, xy=(0, 0.5), xycoords='axes fraction',
                                va='center', ha='left', fontsize=12, color='green')
            axes[i, j].annotate(right, xy=(1, 0.5), xycoords='axes fraction',
                                va='center', ha='right', fontsize=12, color='green')
            axes[i, j].annotate(top, xy=(0.5, 1), xycoords='axes fraction',
                                va='top', ha='center', fontsize=12, color='green')
            axes[i, j].annotate(bottom, xy=(0.5, 0), xycoords='axes fraction',
                                va='bottom', ha='center', fontsize=12, color='green')
    plt.suptitle(f"QC Slices: {fname}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, fname.replace('.nii', '').replace('.gz', '') + '_qc.png')
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
                if fname.endswith('.nii') or fname.endswith('.nii.gz'):
                    nifti_path = os.path.join(modality_dir, fname)
                    img = nib.load(nifti_path)
                    shape = img.shape
                    n_vols = shape[3] if (modality in ['dwi', 'func'] and len(shape) > 3) else 1
                    qc_img_path = plot_nifti_slices(nifti_path, out_dir, n_slices)
                    report_entries.append({
                        'filename': fname,
                        'modality': modality,
                        'dimensions': shape,
                        'n_volumes': n_vols,
                        'qc_img_path': os.path.basename(qc_img_path)
                    })
    # Sort entries by modality and filename
    report_entries.sort(key=lambda x: (x['modality'], x['filename']))
    return report_entries

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
    report_fname = f"sub-{subject_id}_ses-{session_id}_qc_report.html"
    html_path = os.path.join(out_dir, report_fname)
    report_title = f"NIfTI QC Report for {subject_id} {session_id}"
    with open(html_path, "w") as f:
        f.write(f"<html><head><title>{report_title}</title>\n")
        f.write("""
        <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .qc-entry { margin-bottom: 40px; }
        .qc-info { font-size: 1.1em; margin-bottom: 8px; }
        .qc-img { width: 100%; max-width: 1200px; border: 1px solid #ccc; }
        </style>
        """)
        f.write("</head><body>\n")
        f.write(f"<h1>{report_title}</h1>\n")
        for entry in report_entries:
            # Extract subject and session from filename
            parts = entry['filename'].split('_')
            subj = parts[0].replace('sub-', '') if len(parts) > 0 else "unknown"
            sess = parts[1].replace('ses-', '') if len(parts) > 1 else "unknown"
            f.write(f"<div class='qc-entry' data-subject='{subj}' data-session='{sess}' data-modality='{entry['modality']}'>\n")
            f.write(
                f"<div class='qc-info'><b>File Name:</b> {entry['filename']} &nbsp; "
                f"<b>Modality:</b> {entry['modality']} &nbsp; "
                f"<b>Dimensions:</b> {entry['dimensions']} &nbsp; "
                f"<b># Volumes:</b> {entry['n_volumes']}</div>\n"
            )
            f.write(f"<img class='qc-img' src='{entry['qc_img_path']}' alt='{entry['filename']}'>\n")
            f.write("</div>\n")
        # Build dropdown menu for navigation
        subjects = sorted(set(
            entry['filename'].split('_')[0].replace('sub-', '') for entry in report_entries
        ))
        sessions = sorted(set(
            entry['filename'].split('_')[1].replace('ses-', '') for entry in report_entries
        ))
        modalities = sorted(set(entry['modality'] for entry in report_entries))
        # Create dropdown HTML
        f.write("""
        <script>
        function filterQC() {
            var subj = document.getElementById('subjectDropdown').value;
            var sess = document.getElementById('sessionDropdown').value;
            var mod = document.getElementById('modalityDropdown').value;
            var entries = document.getElementsByClassName('qc-entry');
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
        <div id='qc-dropdown-bar' style='position:fixed;top:0;left:0;width:100%;background:#f9f9f9;border-bottom:1px solid #ccc;z-index:1000;padding:12px 0;'>
            <label style='margin-right:20px;'>Subject:
            <select id='subjectDropdown' onchange='filterQC()'>
                <option value='all'>All</option>
        """)
        for s in subjects:
            f.write(f"<option value='{s}'>{s}</option>\n")
        f.write("""
            </select>
            </label>
            <label style='margin-right:20px;'>Session:
            <select id='sessionDropdown' onchange='filterQC()'>
                <option value='all'>All</option>
        """)
        for s in sessions:
            f.write(f"<option value='{s}'>{s}</option>\n")
        f.write("""
            </select>
            </label>
            <label style='margin-right:20px;'>Modality:
            <select id='modalityDropdown' onchange='filterQC()'>
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
        f.write("</body></html>\n")
    print(f"HTML report written to {html_path}")

def main():
    parser = argparse.ArgumentParser(description="Create QC plots and HTML report for NIfTI files of a subject.")
    parser.add_argument("subject_dir", help="Path to sub-<subject_id>/ses-<session_id>/ directory")
    parser.add_argument("out_dir", help="Directory to save QC plots and HTML report")
    parser.add_argument("--n_slices", type=int, default=10, help="Number of slices per orientation")
    args = parser.parse_args()
    report_entries = process_subject(args.subject_dir, args.out_dir, args.n_slices)
    write_html_report(report_entries, args.out_dir)

if __name__ == "__main__":
    main()