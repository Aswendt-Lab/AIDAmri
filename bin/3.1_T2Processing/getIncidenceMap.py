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
    """
    Uses only the last two path elements of inputLocation.
    Replaces slashes with underscores.
    Example:
      /aida/Data/Division/proc/Stroke/ses-PT3
      -> heatmap_Stroke_ses-PT3
    """
    # Normalize path
    input_abs = os.path.abspath(os.path.normpath(inputLocation))

    # Split path into components
    parts = input_abs.split(os.sep)

    # Remove empty elements (important if path starts with /)
    parts = [p for p in parts if p]

    if len(parts) >= 2:
        rel = "_".join(parts[-2:])
    else:
        rel = parts[-1]

    return f"{prefix}{rel}"

def heatMap(incidenceMap, araVol, outputLocation, prefix):
    maxV = int(np.max(incidenceMap))
    fig, axes = plt.subplots(nrows=3, ncols=4)
    t = 1
    for ax in axes.flat:
        im = ax.imshow(np.transpose(np.round(incidenceMap[:, :, t * 16])), cmap='gnuplot', vmin=0, vmax=maxV)
        ax.imshow(np.transpose(araVol[:, :, t * 16]), alpha=0.55, cmap='gray')
        ax.axis('off')
        t = t + 1

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


def incidenceMap2(path_listInc, araTemplate, inputLocation, outputLocation, prefix):
    araDataTemplate = nii.load(araTemplate)
    realAraImg = np.asanyarray(araDataTemplate.dataobj)
    overlaidIncidences = np.zeros_like(realAraImg)
    bar = progressbar.ProgressBar()
    for fileIndex in bar(range(len(path_listInc))):
        dataMRI = nii.load(path_listInc[fileIndex])
        volumeMRI = np.asanyarray(dataMRI.dataobj)

        # Adjusting the volumeMRI data
        volumeMRI[volumeMRI <= 0] = 0
        volumeMRI[volumeMRI > 0] = 1

        overlaidIncidences += volumeMRI

    overlayNII = nii.Nifti1Image(overlaidIncidences, araDataTemplate.affine)
    # remove "heatmap_" from prefix
    name_part = prefix.replace("heatmap_", "", 1)

    output_file = os.path.join(outputLocation, f"incMap_{name_part}.nii.gz")
    nii.save(overlayNII, output_file)
    heatMap(incidenceMap=overlaidIncidences, araVol=realAraImg, outputLocation=outputLocation, prefix=prefix)
    max_overlap = int(np.max(overlaidIncidences))
    print("Maximum number of subjects overlapping at any voxel in the incidence volume:", max_overlap)


def findIncData(path):
    regMR_list = []
    for filename in glob.iglob(os.path.join(path,"*","*",'anat', '*IncidenceData_mask.nii.gz')):
        regMR_list.append(filename)
    return regMR_list


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate an Incidence Map')
    parser.add_argument('-i', '--inputLocation', help='Directory: Brain extracted input data, e.g proc_data folder', required=True)
    parser.add_argument('-o', '--outputLocation', help='Directory: Output location for the heat map', default=None)
    parser.add_argument('-a', '--allenBrainTemplate', help='File: Annotations of Allen Brain', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir, 'lib', 'average_template_50.nii.gz')))

    args = parser.parse_args()

    inputLocation = args.inputLocation
    outputLocation = args.outputLocation
    allenBrainTemplate = args.allenBrainTemplate

    # If no output location is provided → use input directory
    if outputLocation is None:
        outputLocation = inputLocation

    prefix = build_output_prefix(inputLocation)

    if not os.path.exists(inputLocation):
        sys.exit("Error: '%s' is not an existing directory." % (inputLocation,))

    os.makedirs(outputLocation, exist_ok=True)

    if not os.path.exists(allenBrainTemplate):
        sys.exit("Error: '%s' is not an existing file." % (allenBrainTemplate,))

    regInc_list = findIncData(inputLocation)

    if len(regInc_list) < 1:
        sys.exit("Error: No masked strokes found in the provided directory.")

    print("'%i' folders are part of the incidence map." % (len(regInc_list),))
    incidenceMap2(regInc_list, allenBrainTemplate, inputLocation, outputLocation, prefix)
    sys.exit(0)
