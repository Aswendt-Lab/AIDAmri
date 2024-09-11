The Atlas-based Imaging Data Analysis Pipeline for structural and functional MRI of the rodent brain (AIDAmri) represents an automated tool specifically designed for standardized data handling and pre-processing [1].
## Details of the Processing Steps Applied to the Multiverse Dataset

### 1. Data Inspection and Manual Re-orientation
As the data was provided in BIDS format, no conversion to NIfTI was necessary. However, some manual re-orientation was required to ensure the images met the input criteria for AIDA<em>mri</em>. Specifically, the data had to be adjusted to the RAS (Neurological convention) orientation, ensuring correct alignment and flipping to ensure compatibility with the pipeline.

### 2. Registration

The registration is performed in two distinct steps to ensure precise alignment. 

2.1. **Affine Transformation**: 
   - The initial step involves an affine transformation applied specifically to T2-weighted (T2) MRI images. This transformation accounts for basic geometric adjustments such as scaling, rotation, translation, and shearing. The affine registration is performed using NiftyReg's symmetric block-matching approach, which ensures a global alignment of the MRI data with the SIGMA rat brain MRI template. The template serves as an intermediate reference that closely correlates with the anatomical features captured in T2-weighted images.

2.2. **Non-linear Transformation**:
   - After the affine registration, a non-linear transformation is applied to address more complex deformations that cannot be captured by affine transformations alone. This step is crucial for accurately mapping subcortical brain structures, particularly in cases involving significant brain deformation, such as stroke-induced lesions. The non-linear transformation is performed on the MRI template, which is already aligned with the T2-weighted images. The resulting transformation is then applied to align the data with the ARA.

2.3. **Transfer of Transformations to fMRI**: 
   - Once the affine and non-linear transformations are determined for the T2-weighted images, these transformations are subsequently applied to the functional MRI datasets. This approach leverages the high anatomical contrast in T2-weighted images to ensure accurate registration for more noisy and lower resolution data.

### 3. fMRI Processing

3.1. **Motion Correction**:
   - fMRI data undergoes slice-wise motion correction to address any displacements that occur during scanning.

3.2. **Smoothing and Filtering**:
   - A spatial filter is applied to smooth the data, focusing on the x-y plane to account for anisotropy in voxel sizes. A high-pass filter with a cut-off frequency of 0.01 Hz is used to remove low-frequency noise, improving the signal quality.

3.4. **Time-Series Extraction**:
   - For each region of interest as defined by the ARA, the mean intensity of the voxels is calculated over time, producing a time series that reflects the region's activity.

3.5. **Functional Connectivity Analysis**:
   - Correlation of BOLD signals between brain regions is calculated, providing insights into the functional connectivity across the brain. This analysis helps identify networks of interacting regions based on their synchronized activity patterns.

### Multiverse-specific output

### Tools used by AIDAmri

| Tool          | Version         | Description                                              | 
|---------------|-----------------|----------------------------------------------------------|
| [**Docker**](https://docs.docker.com/engine/install/)                                   | latest                      | Software to execute the AIDAmri Docker file
| [**Brkraw**](https://github.com/brkraw/bruker@e27c5039c9c3a84ce7cd19c9627360e5a31b4ebc) | specific commit | Software to convert Bruker into BIDS            
| [**NiftyReg**](https://sourceforge.net/projects/niftyreg/)                              | 1.5.55                       | Software library for medical image registration     
| [**FSL**](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL)                                   | 5.0.11                      | Library of analysis tools for brain imaging data  
| [**DSI Studio**](http://dsi-studio.labsolver.org/)                                      | 2023.07.08                               | Diffusion MRI analysis tool                                                  

### Reference
[1]: Processing pipeline for Atlas-based Imaging Data Analysis (AIDA) of structural and functional mouse brain MRI. . Pallast N, Diedenhofen M, Blaschke S, Wieters F, Wiedermann D, Hoehn M,  Fink GR, Aswendt M. Front Neuroinform . 2019 Jun 4;13:42. doi: 10.3389/fninf.2019.00042. 
