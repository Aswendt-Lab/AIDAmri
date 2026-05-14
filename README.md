![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/Aswendt-Lab/AIDAmri/docker-image.yml) ![Static Badge](https://img.shields.io/badge/Docker_image-11.97_GB-blue) [![Static Badge](https://img.shields.io/badge/data_structure-BIDS-yellow)](https://bids.neuroimaging.io/news.html) [![Static Badge](https://img.shields.io/badge/Niftyreg-CBSI-orange)](https://github.com/KCL-BMEIS/niftyreg) [![Static Badge](https://img.shields.io/badge/DSI--Studio-2023-orange)](https://dsi-studio.labsolver.org/) [![Static Badge](https://img.shields.io/badge/FSL-5.0.11-orange)]([https://dsi-studio.labsolver.org/](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)) ![Static Badge](https://img.shields.io/badge/Python-3.7-orange)

[1.2]: http://i.imgur.com/wWzX9uB.png
[1]: http://www.twitter.com/AswendtMarkus
<!--social icon from https://github.com/carlsednaoui/gitsocial -->

<img align="left" src="https://github.com/maswendt/AIDAmri/blob/master/AIDA_Logo.png" width="120">
<h1>AIDA<i>mri</i></h1>

Atlas-based Imaging Data Analysis Pipeline (AIDA) for structural and functional MRI of the mouse brain
<br/>
## Key Features of AIDA<em>mri</em>

1. **Automated Preprocessing**  
   It performs tasks like image re-orientation, bias-field correction, and brain extraction with minimal user input required.

2. **Atlas-Based Registration**  
   AIDA<em>mri</em> uses the **Allen Mouse Brain Reference Atlas** for accurate region-based analysis of MRI data, allowing researchers to compare results across different studies efficiently. A modified atlas version with larger labels to better match MRI resolution is provided. Users can define specific **regions of interest (ROIs)** for analysis, such as stroke lesions.

3. **Modular Design**  
   The pipeline is developed in Python, making it cross-platform and open-source, allowing for easy integration and modification.

4. **Validation**  
   The pipeline was validated with different MRI datasets, including those involving stroke models, demonstrating its robustness even in the presence of significant brain deformations.

5. **Functional and Structural Connectivity Analysis**  
   The output of the pipeline includes connectivity matrices that can be used for further analysis of brain network changes in health and disease.

<p align="center">
  <img src="https://github.com/maswendt/AIDAmri/blob/master/AIDAmri_drawing.png" style="max-width: 100%; height: auto;">
</p>

Pipeline overview from [Pallast et al.](https://doi.org/10.3389/fninf.2019.00042)



## Version history

[Information latest Version 2.0](https://github.com/maswendt/AIDAmri/releases/tag/v2.0)

[**Manual**](https://github.com/maswendt/AIDA/blob/master/manual.pdf)

[Information about Version 1.2 (Docker stable release)](https://github.com/maswendt/AIDAmri/releases/tag/v1.2)
<br/>
[Information about Version 1.1.1 (Docker pre-release)](https://github.com/maswendt/AIDAmri/releases/tag/1.1.1)
<br/>
[Information about Version 1.1 (Stable)](https://github.com/maswendt/AIDAmri/releases/tag/v1.1)
<br/>
[Information about Version 1.0](https://github.com/maswendt/AIDAmri/releases/tag/v1.0)

<h3><b>Important note: read this before you install AIDAmri for the first time</h3></b>

We fully moved to the containerized version of AIDAmri via [Docker](https://docs.docker.com/get-docker/). All information can be found in the manual above. Please report issues and bugs directly in the issue section of this repository or at gitter (Link below in the contact section).

## BRANCHES

AIDAmri is organized into multiple branches to support development, collaboration, and species-specific adaptations:

- **`main`** – the stable branch containing officially released and validated versions of AIDAmri for mice.  
- **`open-dev`** – the public development branch that can be used by external contributors to implement code modifications, enhancements, or bug fixes.  
  *Researchers and developers are welcome to fork the repository, work within the `open-dev` branch, and submit pull requests for review.*  
- **`rat`** – a dedicated branch for rat MRI data, including modified atlas and template versions optimized for rat brain imaging and analysis.  

Each branch is continuously synchronized to ensure compatibility with the core AIDAmri framework and Docker-based environment.

## EXAMPLE FILES

Download [**here**](https://gin.g-node.org/Aswendt_Lab/testdata_AIDA) (you probably have to clone the dataset from the gin repo. The files are annexed files, also use the raw_data folder as the test data).\
Mouse MRI data, acquired with Bruker 9.4T - cryo coil setup: adult C57BL7/6 mouse, 
T2-weighted (anatomical scan),
DTI (structural connectivity scan),
rs-fMRI (functional connectivity scan).

## Data Format and Orientation Requirements

AIDAmri supports data processing exclusively for datasets in NIfTI (.nii/.nii.gz) or Bruker formats. To ensure accurate registration and reproducible results, 
all input data for preprocessing must be in <ins>**LIP (Left-Inferior-Posterior)**</ins> orientation.
Furthermore, the image header information must be consistent with the physical orientation of the data array. 
Any mismatch between the header orientation and the actual voxel layout can lead to registration errors or incorrect alignment with the atlas. It is therefore strongly recommended to verify and, if necessary, correct the header orientation. 
Please use FSL eyes for visual inspection and fslhd for checking the header information. FSL is already installed inside AIDAmri. More Information about FSL can be find [here](https://fsl.fmrib.ox.ac.uk/fsl/docs/).
If your data is in a different orientation than LIP please use our ReorientBatch.py script in the helpertools folder. The script should be used after convert2Nifti script and can reorient the whole proc_data folder. 
It is important that the folder contains only the NIFTI files to be reoriented. The folder must not contain any NIFTI files that have already been processed. 
Furthermore, please note that after reorientation, tools such as Fiji or other tools that do not read the header of a NIFTI file will display the images only as the data was saved after reorientation. 
For this reason, we recommend FSL Eyes, as this tool provides more information about the orientation. 


## Troubleshooting / Common Issues

<details>
<summary><strong>This section lists frequently encountered problems when using AIDAmri and possible solutions.</strong></summary>
	
If your problem is not listed here, please use our Gitter Chat or open an issue on GitHub and include:
- OS
- Docker version
- Command used
- Full error log


---

<details>
<summary><strong>General debugging tips</strong></summary>

Always process **T2 data** first

Visually inspect outputs after:

- Bias correction

- Brain extraction

- Registration

Check logs instead of relying only on exit codes

</details>

<details>
<summary><strong>NIfTI orientation problems</strong></summary>
	
**Symptom**

- Atlas is flipped or registration fails

- Brain appears flipped or rotated


**Explanation**

AIDAmri expects LIP orientation for preprocessing.

**Solution**

We recommend verifying the NIfTI header information and the actual image orientation using the FSL tools `fsleyes` and `fslhd`.

In cases where the images are not in a consistent LIP orientation, the provided reorientation script should be applied **after** `convert2nifti` and **prior to preprocessing**.:

```
python ReorientBatch.py -i proc_data -o proc_data_reoriented
```

</details>


<details>
<summary><strong>Docker build takes a very long time, finishes after under 1 min or aborts</strong></summary>
	
**Symptom**

docker build hangs for several minutes or fails during apt-get / FSL installation

**Explanation**

- Possible causes

- Slow internet connection

- Docker build cache corrupted

- Insufficient disk space

**Solution**
```
docker system prune -a

docker build --no-cache -t aidamri:latest .
```

</details>

<details>
<summary><strong>A is singular, uses pseudoinverse in updateB</strong></summary>
	
**Symptom**

`Warning: A is singular, uses pseudoinverse in updateB`

**Explanation**

This warning originates from the MICO bias-field correction and usually indicates:

- Low SNR

- Strong intensity inhomogeneities

- Empty or corrupted slices

**Solution**

- Usually safe to ignore if preprocessing finishes successfully

- Check input image for NaNs, zero-only slices or low voxel values

- Verify correct image orientation before preprocessing

</details>

<details>
<summary><strong>Brain extraction fails (SVD did not converge)</strong></summary>
	
**Symptom**

`Error in brain extraction`

`SVD did not converge`

**Explanation**

- Corrupted NIfTI header

- NaN or Inf values in the image

- Extreme bias-field artifacts

**Solution**

- Re-run bias-field correction

- Reorient the image manually before preprocessing

- Check header consistency:
```
fslhd input.nii.gz
```

</details>

</details>


## ARA CREATOR
[Matlab script](https://github.com/maswendt/AIDAmri/ARA) to generate a custom version of the Allen Mouse Brain Atlas.

[<h3><b>CONTACT</h3></b>]
If you encounter problems, report directly in [![Gitter](https://badges.gitter.im/AIDA_tools/community.svg)](https://gitter.im/AIDA_tools/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

or 

join our Open Office Hour - each Thursday 3:00 pm (UTC+2) [![Zoom](https://img.shields.io/badge/Zoom-2D8CFF?style=for-the-badge&logo=zoom&logoColor=white)](https://uni-frankfurt.zoom-x.de/j/63112745009?pwd=JBTjMVbuaTw9cZvFnppTwCPjGdQEyx.1)


For all other inquiries: Markus Aswendt (aswendtATmed.uni-frankfurt.de)

<h3><b>LICENSE/CITATION</h3></b>
GNU General Public License v3.0
<br/>
<br/>
If you use our software or modify parts of it and use it in other ways, please cite: 
<br/>
<br/>

*Pallast N, Diedenhofen M, Blaschke S, Wieters F, Wiedermann D, Hoehn M, Fink GR, Aswendt M. Processing Pipeline for Atlas-Based Imaging Data Analysis of Structural and Functional Mouse Brain MRI (AIDAmri). Front Neuroinform. 2019 Jun 4;13:42.[doi: 10.3389/fninf.2019.00042.](https://doi.org/10.3389/fninf.2019.00042)*
___
<details>
<summary>REFERENCES</summary></b>

+ Brain Connectivity Toolbox
    + [M. Rubinov and O. Sporns (2010). Complex Network Measures of Brain Connectivity: Uses 
and Interpretations. NeuroImage 52 (3), 1059–69.](https://www.sciencedirect.com/science/article/abs/pii/S105381190901074X)
+ Allen Mouse Brain Reference Atlas
    + [Wang et al. (2020). The Allen Mouse Brain Common Coordinate Framework: A 3D Reference Atlas. Cell 181 (4), 936-953.](https://pubmed.ncbi.nlm.nih.gov/32386544/)
+ Niftyreg
    + [Ourselin, et al. (2001). Reconstructing a 3D structure from serial
histological sections. Image and Vision Computing, 19(1-2), 25–31.](https://www.sciencedirect.com/science/article/pii/S0262885600000524)
    + [Modat, et al. (2014). Global image registration using a symmetric block-
matching approach. Journal of Medical Imaging, 1(2), 024003–024003.](https://www.ncbi.nlm.nih.gov/pubmed/26158035)
    + [Rueckert, et al.. (1999). Nonrigid registration using free-form
deformations: Application to breast MR images. IEEE Transactions on Medical
Imaging, 18(8), 712–721.](https://ieeexplore.ieee.org/document/796284)
    + [Modat, et al. (2010). Fast free-form deformation using graphics processing
units. Computer Methods And Programs In Biomedicine,98(3), 278–284.](https://www.ncbi.nlm.nih.gov/pubmed/19818524)
+ FSL
    + [M.W. Woolrich, S. Jbabdi, B. Patenaude, M. Chappell, S. Makni, T. Behrens, C. Beckmann, M. Jenkinson, S.M. Smith. Bayesian analysis of neuroimaging data in FSL. NeuroImage, 45:S173-86, 2009](https://www.ncbi.nlm.nih.gov/pubmed/19059349)
    + [S.M. Smith, M. Jenkinson, M.W. Woolrich, C.F. Beckmann, T.E.J. Behrens, H. Johansen-Berg, P.R. Bannister, M. De Luca, I. Drobnjak, D.E. Flitney, R. Niazy, J. Saunders, J. Vickers, Y. Zhang, N. De Stefano, J.M. Brady, and P.M. Matthews. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 23(S1):208-19, 2004](https://www.sciencedirect.com/science/article/pii/S1053811904003933?via%3Dihub)
    + [M. Jenkinson, C.F. Beckmann, T.E. Behrens, M.W. Woolrich, S.M. Smith. FSL. NeuroImage, 62:782-90, 2012](https://www.sciencedirect.com/science/article/pii/S1053811911010603?via%3Dihub) 
+ DSIstudio
    + [Yeh, Fang-Cheng, et al. Deterministic diffusion fiber tracking improved by quantitative anisotropy. (2013): e80713. PLoS ONE 8(11)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080713)
</details>
