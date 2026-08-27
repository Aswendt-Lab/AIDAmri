import os
import sys
import nibabel as nii
import glob
import numpy as np
import progressbar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Fonts & Text Display ---
matplotlib.rcParams['svg.fonttype'] = 'none'     #text remains editable in SVG
matplotlib.rcParams['pdf.fonttype'] = 42         # Editable text in PDF (Type 42)

def build_output_prefix(inputLocation, prefix="heatmap_"):
    input_abs = os.path.abspath(os.path.normpath(inputLocation))
    input_name = os.path.basename(input_abs)
    if not input_name:
        sys.exit("Error: Input location must not be the filesystem root.")
    return f"{prefix}{input_name}"


def validate_heatmap_name(heatmap_name):
    heatmap_name = heatmap_name.strip()
    if not heatmap_name:
        sys.exit("Error: Heatmap name must not be empty.")
    if os.path.basename(heatmap_name) != heatmap_name:
        sys.exit("Error: Heatmap name must be a file name, not a path.")
    return heatmap_name


def validate_session(session):
    session = session.strip()
    if not session:
        sys.exit("Error: Session must not be empty.")
    if os.path.basename(session) != session:
        sys.exit("Error: Session must be a folder name, not a path.")
    return session


def heatMap(incidenceMap, araVol, outputLocation, prefix):
    maxV = int(np.max(incidenceMap))
    fig, axes = plt.subplots(nrows=3, ncols=4)
    z_slices = np.linspace(0, incidenceMap.shape[2] - 1, 12, dtype=int)
    for ax, z_index in zip(axes.flat, z_slices):
        im = ax.imshow(np.transpose(np.round(incidenceMap[:, :, z_index])), cmap='gnuplot', vmin=0, vmax=maxV)
        ax.imshow(np.transpose(araVol[:, :, z_index]), alpha=0.55, cmap='gray')
        ax.axis('off')

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    bounds = np.linspace(0, maxV, maxV + 1)
    cbar = fig.colorbar(im, cax=cbar_ax, format='%1i', ticks=bounds)
    cbar.ax.tick_params(labelsize=14)

    # Save the heatmap instead of showing
    output_file = os.path.join(outputLocation, f"{prefix}.png")
    fig.savefig(output_file, dpi=300, bbox_inches="tight")

    # Save heatmap as PDF
    output_pdf = os.path.join(outputLocation, f"{prefix}.pdf")
    fig.savefig(output_pdf, bbox_inches="tight")

    # Save heatmap as SVG (vector graphics)
    output_svg = os.path.join(outputLocation, f"{prefix}.svg")
    fig.savefig(output_svg, bbox_inches="tight")

    plt.close(fig)


def incidenceMap2(path_listInc, araTemplate, outputLocation, prefix):
    araDataTemplate = nii.load(araTemplate)
    realAraImg = np.asanyarray(araDataTemplate.dataobj)
    overlaidIncidences = np.zeros(realAraImg.shape, dtype=np.uint16)
    bar = progressbar.ProgressBar()
    for fileIndex in bar(range(len(path_listInc))):
        dataMRI = nii.load(path_listInc[fileIndex])
        volumeMRI = np.asanyarray(dataMRI.dataobj).copy()

        if volumeMRI.shape != overlaidIncidences.shape:
            sys.exit(
                "Error: Shape mismatch for '%s'. Mask shape is %s, template shape is %s." %
                (path_listInc[fileIndex], volumeMRI.shape, overlaidIncidences.shape,)
            )
        if not np.allclose(dataMRI.affine, araDataTemplate.affine, atol=1e-4):
            sys.exit(
                "Error: Affine mismatch for '%s'. Mask and template are not in the same grid." %
                (path_listInc[fileIndex],)
            )

        # Adjusting the volumeMRI data
        volumeMRI[volumeMRI <= 0] = 0
        volumeMRI[volumeMRI > 0] = 1
        volumeMRI = volumeMRI.astype(overlaidIncidences.dtype)

        overlaidIncidences += volumeMRI

    overlayNII = nii.Nifti1Image(overlaidIncidences, araDataTemplate.affine)
    # remove "heatmap_" from prefix
    name_part = prefix.replace("heatmap_", "", 1)

    output_file = os.path.join(outputLocation, f"incMap_{name_part}.nii.gz")
    nii.save(overlayNII, output_file)
    heatMap(incidenceMap=overlaidIncidences, araVol=realAraImg, outputLocation=outputLocation, prefix=prefix)
    max_overlap = int(np.max(overlaidIncidences))
    print("Maximum number of subjects overlapping at any voxel in the incidence volume:", max_overlap)


def findIncData(path, session):
    search_pattern = os.path.join(
        path,
        "*",
        session,
        "anat",
        "IncidenceData",
        "*IncidenceData_Lesion_mask.nii.gz",
    )
    return sorted(glob.glob(search_pattern))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate an Incidence Map')
    parser.add_argument('-i', '--inputLocation', help='Directory: Brain extracted input data, e.g proc_data folder', required=True)
    parser.add_argument('--session', help='Session folder to include, e.g. ses-PT3', required=True)
    parser.add_argument('-o', '--outputLocation', help='Directory: Output location for the heat map', default=None)
    parser.add_argument('-n', '--heatmapName', help='Optional output name for the heatmap files', default=None)
    parser.add_argument('-a', '--allenBrainTemplate', help='File: Annotations of Allen Brain', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir, 'lib', 'average_template_50.nii.gz')))

    args = parser.parse_args()

    inputLocation = args.inputLocation
    outputLocation = args.outputLocation
    allenBrainTemplate = args.allenBrainTemplate
    heatmapName = args.heatmapName
    session = validate_session(args.session)

    # If no output location is provided → use input directory
    if outputLocation is None:
        outputLocation = inputLocation

    if heatmapName is None:
        prefix = f"{build_output_prefix(inputLocation)}_{session}"
    else:
        prefix = validate_heatmap_name(heatmapName)

    if not os.path.exists(inputLocation):
        sys.exit("Error: '%s' is not an existing directory." % (inputLocation,))

    os.makedirs(outputLocation, exist_ok=True)

    if not os.path.exists(allenBrainTemplate):
        sys.exit("Error: '%s' is not an existing file." % (allenBrainTemplate,))

    regInc_list = findIncData(inputLocation, session)

    if len(regInc_list) < 1:
        sys.exit("Error: No masked strokes found for session '%s' in the provided directory." % (session,))

    print("'%i' folders from session '%s' are part of the incidence map." % (len(regInc_list), session,))
    incidenceMap2(regInc_list, allenBrainTemplate, outputLocation, prefix)
    sys.exit(0)
