import os
import sys
import nibabel as nii
import glob
import numpy as np
import progressbar
import matplotlib
import matplotlib.pyplot as plt

# --- Fonts & Text Display ---
matplotlib.rcParams['svg.fonttype'] = 'none'     #text remains editable in SVG
matplotlib.rcParams['pdf.fonttype'] = 42         # Editable text in PDF (Type 42)

def heatMap(incidenceMap, araVol, outputLocation):
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
    output_file = os.path.join(outputLocation, 'heatMap.png')
    plt.savefig(output_file)

    # Save heatmap as PDF
    output_pdf = os.path.join(outputLocation, 'heatMap.pdf')
    plt.savefig(output_pdf)

    # Save heatmap as SVG (vector graphics)
    output_svg = os.path.join(outputLocation, 'heatMap.svg')
    plt.savefig(output_svg)

    plt.close()


def incidenceMap2(path_listInc, araTemplate, inputFile, outputLocation):
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
    output_file = os.path.join(outputLocation, 'incMap.nii.gz')
    nii.save(overlayNII, output_file)
    heatMap(incidenceMap=overlaidIncidences, araVol=realAraImg, outputLocation=outputLocation)
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
    parser.add_argument('-i', '--inputFile', help='Directory: Brain extracted input data, e.g proc_data folder', required=True)
    parser.add_argument('-o', '--outputLocation', help='Directory: Output location for the heat map', required=True)
    parser.add_argument('-a', '--allenBrainTemplate', help='File: Annotations of Allen Brain', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir, 'lib', 'average_template_50.nii.gz')))

    args = parser.parse_args()

    inputFile = args.inputFile
    outputLocation = args.outputLocation
    allenBrainTemplate = args.allenBrainTemplate

    if not os.path.exists(inputFile):
        sys.exit("Error: '%s' is not an existing directory." % (inputFile,))

    if not os.path.exists(outputLocation):
        sys.exit("Error: '%s' is not an existing directory." % (outputLocation,))

    if not os.path.exists(allenBrainTemplate):
        sys.exit("Error: '%s' is not an existing file." % (allenBrainTemplate,))

    regInc_list = findIncData(inputFile)

    if len(regInc_list) < 1:
        sys.exit("Error: No masked strokes found in the provided directory.")

    print("'%i' folders are part of the incidence map." % (len(regInc_list),))
    incidenceMap2(regInc_list, allenBrainTemplate, inputFile, outputLocation)
    sys.exit(0)
