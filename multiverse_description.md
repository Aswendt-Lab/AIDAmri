The Atlas-based Imaging Data Analysis Pipeline for structural and functional MRI of the rodent brain (AIDAmri) represents an automated tool specifically designed for standardized data handling and pre-processing [1].
## Details of the Processing Steps Applied to the Multiverse Dataset

### Data Inspection and Manual Re-orientation
As the data was provided in BIDS format, no conversion to NIfTI was necessary. However, some manual re-orientation was required to ensure the images met the input criteria for AIDA<em>mri</em>. Specifically, the data had to be adjusted to the RAS (Neurological convention) orientation, ensuring correct alignment and flipping to ensure compatibility with the pipeline.

### Registration

#### Affine and Non-linear Transformations
The registration is performed in two distinct steps to ensure precise alignment:

1. **Affine Transformation**: The initial step involves an affine transformation applied specifically to T2-weighted (T2) MRI images. This transformation accounts for basic geometric adjustments such as scaling, rotation, translation, and shearing. The affine registration is performed using NiftyReg's symmetric block-matching approach, which ensures a global alignment of the MRI data with a custom-developed MRI template (MTPL). The MTPL serves as an intermediate reference that closely correlates with the anatomical features captured in T2-weighted images.

2. **Non-linear Transformation**: After the affine registration, a non-linear transformation is applied to address more complex deformations that cannot be captured by affine transformations alone. This step is crucial for accurately mapping subcortical brain structures, particularly in cases involving significant brain deformation, such as stroke-induced lesions. The non-linear transformation is performed on the MRI template (MTPL), which is already aligned with the T2-weighted images. The resultant transformation is then applied to align the data with the ARA.

#### Intermediate MRI Template (MTPL)
The use of an MRI template (MTPL) as an intermediate step is a key feature of AIDA<em>mri</em>’s registration process. The MTPL is constructed from multiple T2-weighted MRI datasets of healthy C57BL6 mice, providing a high-resolution average image that reflects the common anatomical features across different mice. This template is crucial because it bridges the gap between individual mouse MRI data and the ARA, which was originally created from two-photon microscopy images and has lower correlation with MRI data.

#### Transfer of Transformations to DTI and fMRI
Once the affine and non-linear transformations are determined for the T2-weighted images, these transformations are subsequently applied to the diffusion tensor imaging (DTI) and resting-state functional MRI (rs-fMRI) datasets. This approach leverages the high anatomical contrast in T2-weighted images to ensure accurate registration for other imaging modalities, which may have lower contrast or be more affected by artifacts. By transferring the transformations derived from T2-weighted images, AIDA<em>mri</em> ensures consistent and accurate alignment across different types of MRI data, facilitating comprehensive structural and functional analysis.

### rs-fMRI Processing

1. **Physiological Recording Correction**:
   - **Method**: Respiratory signals recorded during the scan are used to correct for artifacts related to breathing. The code identifies and adjusts for these fluctuations, ensuring accurate signal analysis.

2. **Motion Correction**:
   - **Procedure**: Similar to the DTI motion correction, rs-fMRI data undergoes slice-wise motion correction to address any displacements that occur during scanning.

3. **Smoothing and Filtering**:
   - **Method**: A spatial filter is applied to smooth the data, focusing on the x-y plane to account for anisotropy in voxel sizes. A high-pass filter with a cut-off frequency of 0.01 Hz is used to remove low-frequency noise, improving the signal quality.

#### 3.3 fMRI Activity

1. **Time-Series Extraction**:
   - **Procedure**: For each region of interest as defined by the ARA, the mean intensity of the voxels is calculated over time, producing a time series that reflects the region's activity.

2. **Functional Connectivity Analysis**:
   - **Method**: Correlation of BOLD signals between brain regions is calculated, providing insights into the functional connectivity across the brain. This analysis helps identify networks of interacting regions based on their synchronized activity patterns.

[1]: Processing pipeline for Atlas-based Imaging Data Analysis (AIDA) of structural and functional mouse brain MRI. . Pallast N, Diedenhofen M, Blaschke S, Wieters F, Wiedermann D, Hoehn M,  Fink GR, Aswendt M. Front Neuroinform . 2019 Jun 4;13:42. doi: 10.3389/fninf.2019.00042. 
