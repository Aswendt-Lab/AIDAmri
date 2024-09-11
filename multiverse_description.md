The Atlas-based Imaging Data Analysis Pipeline for structural and functional MRI of the rodent brain (AIDAmri) represents an automated tool specifically designed for standardized data handling and pre-processing [1].
## Details of the Processing Steps Applied to the Multiverse Dataset

### 1. Data Inspection and Manual Re-orientation
As the data was provided in BIDS format, no conversion to NIfTI was necessary. However, some manual re-orientation was required to ensure the images met the input criteria for AIDA<em>mri</em>. Specifically, the data had to be adjusted to the RAS (Neurological convention) orientation, ensuring correct alignment and flipping to ensure compatibility with the pipeline.

### 2. Registration

2.1. **Affine Transformation**: 
   - The initial step involves an affine transformation applied specifically to T2-weighted (T2) MRI images. This transformation accounts for basic geometric adjustments such as scaling, rotation, translation, and shearing. The affine registration is performed using NiftyReg's symmetric block-matching approach, which ensures a global alignment of the MRI data with the SIGMA rat brain MRI template. The template serves as an intermediate reference that closely correlates with the anatomical features captured in T2-weighted images.

2.2. **Non-linear Transformation**:
   - After the affine registration, a non-linear transformation is applied to address more complex deformations that cannot be captured by affine transformations alone. This step is crucial for accurately mapping subcortical brain structures, particularly in cases involving significant brain deformation, such as stroke-induced lesions. The non-linear transformation is performed on the MRI template, which is already aligned with the T2-weighted images. The resulting transformation is then applied to align the data with the ARA.

2.3. **Transfer of Transformations to fMRI**: 
   - Once the affine and non-linear transformations are determined for the T2-weighted images, these transformations are subsequently applied to the functional MRI datasets. This approach leverages the high anatomical contrast in T2-weighted images to ensure accurate registration for more noisy and lower resolution data.

### 3. fMRI Processing

3.1 **Motion Correction, Smoothing and Filterin**
   - The fMRI data undergoes slice-wise motion correction. This correction reduces artifacts caused by head movements during scanning, ensuring that each slice is properly aligned over time.  To reduce noise, a spatial filter is applied, focusing on the x-y plane to address anisotropy in voxel sizes. A high-pass filter with a cut-off frequency of 0.01 Hz is also used to eliminate low-frequency noise, thereby improving the signal quality. The script applies a series of steps to mask non-brain areas, smooth the data using SUSAN (a non-linear noise reduction algorithm), and apply thresholding to remove low-intensity voxels. 

3.3 **Brain Extraction (BET) and Bias Field Correction**
   - FSL bet is used to isolate brain tissue from non-brain tissue. ANTs N4BiasFieldCorrection corrects intensity inhomogeneities, ensuring uniform signal intensity across the brain for better registration results. 

3.4. **Time-Series Extraction**:
   - For each region of interest as defined by the ARA, the mean intensity of the voxels is calculated over time, producing a time series that reflects the region's activity.

3.5. **Functional Connectivity Analysis**:
   - Correlation of BOLD signals between brain regions is calculated, providing insights into the functional connectivity across the brain. This analysis helps identify networks of interacting regions based on their synchronized activity patterns.

### 4. Graphical Summary
``` mermaid
flowchart TD
    A1["Anatomical T2(x)"] --> B1["Pre-processing"]
    B1 --> D1["Registration"]
  
    A3["Functional fMRI(x)"] --> B3["Pre-processing"]
    B3 --> C3["Regression"]
    C3 --> D3["Registration"]
    D3 --> E3["SIGMA region extraction"]
    
    E3 --> F8["Functional timeseries data"]

    MRI["Multiverse dataset"] --> A1
    MRI --> A3

    D1 -.-> D3



```

### 5. Multiverse-specific output
5.1. **Fully preprocessed EPI in SIGMA space ready for seed-based analysis**
5.2 **The motion parameters**
5.3 **A single frame or temporal mean EPI in SIGMA space with original contrast**
5.4 **Shell script used or a detailed description of each preprocessing step. Including the version for each tool that you used.

5.5 **A list of people to add to the author list (as many as you want), with affiliation, orcid, email, and conflict of interests.**
No conflict of interest. 
Aref Kalantari, Daniel Mertens, and Markus Aswendt. University of Cologne, Faculty of Medicine and University Hospital Cologne, Department of Neurology, Cologne, Germany



### Tools used by AIDAmri

| Tool          | Version         | Description                                              | 
|---------------|-----------------|----------------------------------------------------------|
| [**Docker**](https://docs.docker.com/engine/install/)                                   | latest                      | Software to execute the AIDAmri Docker file
| [**Brkraw**](https://github.com/brkraw/bruker@e27c5039c9c3a84ce7cd19c9627360e5a31b4ebc) | specific commit | Software to convert Bruker into BIDS            
| [**NiftyReg**](https://sourceforge.net/projects/niftyreg/)                              | 1.5.55                       | Software library for medical image registration     
| [**ANTS**](https://github.com/ANTsX/ANTs)                                               | latest                      | N4BiasFieldCorrection algorithm
| [**FSL**](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL)                                   | 5.0.11                      | Library of analysis tools for brain imaging data  
| [**DSI Studio**](http://dsi-studio.labsolver.org/)                                      | 2023.07.08                               | Diffusion MRI analysis tool     


### Reference
[1]: Processing pipeline for Atlas-based Imaging Data Analysis (AIDA) of structural and functional mouse brain MRI. . Pallast N, Diedenhofen M, Blaschke S, Wieters F, Wiedermann D, Hoehn M,  Fink GR, Aswendt M. Front Neuroinform . 2019 Jun 4;13:42. doi: 10.3389/fninf.2019.00042. 
