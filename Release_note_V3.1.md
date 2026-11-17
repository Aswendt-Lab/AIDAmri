# AIDAmri – Changes since v3.0

This document summarizes changes in the v3.1 since the published [v3.0 release](https://github.com/Aswendt-Lab/AIDAmri/releases/tag/v3.0). Features already present in v3.0 are not listed as new additions.

## Highlights

- Added MRtrix3 MP-PCA DWI denoising and improved Patch2Self outputs, including residual images.
- Replaced minimum projection with median projection for DWI/fMRI reference images and made spatial smoothing optional.
- Revised DTI atlas registration using composed transformations.
- Added stroke exclusion, batch processing and cohort CSVs for DTI ROI extraction.
- Expanded DSI connectivity outputs, MATLAB compatibility and cohort extraction.
- Added per-subject and per-session batch parameter overrides through an exemption CSV.
- Expanded quality reports, NIfTI overlays and group filters.
- Added output tracking and targeted resets of processing stages.
- Expanded Git provenance, logging and test coverage.
- Added cohort aggregation and plots of affected brain regions.

## Docker and Dependencies

- Docker now downloads the Ubuntu 22.04 CPU build of DSI Studio `2025.04.16`; library and executable paths use `dsi-studio-cpu`. Build arguments configure the download URL and archive name.
- Added MRtrix3 `3.0.3-1` installation and a build-time `dwidenoise` check.
- Bundled `bet4animal` from `bin/common/bet4animal`, made it executable and removed the additional build-time download. Added `bin/common` to the container PATH and updated the shared BET invocation.
- Added a build stage that generates `/aida/build/AIDAmri_git_information.txt`. The final stage receives the provenance text file without the Git directory.
- Updated `.dockerignore` to include `.git` in the build context for metadata capture.
- Added `install/write_aidamri_git_info.sh` to capture commit, branch, working tree status, author, Git user and repository/remote information. The host Git user can be passed through `AIDAMRI_GIT_CONFIG_USER`.
- Added `install/aidamri_entrypoint.sh` to copy build information to `/aida/DATA` at container startup. The Dockerfile creates the required directories and sets the entrypoint and default Bash command.
- Changed the Seaborn requirement from `==0.11.2` to `>=0.13.2` and its constraint to `0.13.2`.
- Updated the Pillow constraint from `12.2.0` to `12.3.0`.
- Suppressed FSL installer progress displays when no terminal is available or `AIDAMRI_DISABLE_SPINNER=1` is set.
- Ubuntu 22.04, Python 3.10 and `constraints.txt` were already present in v3.0. Some release note additions only clarify these existing details.

## Preprocessing

- Integrated shared output tracking and script logging into T2, DWI, fMRI and T2map preprocessing.
- Original NIfTIs are moved to `brkraw`, then copied back to the input path as working copies.
- Renamed the CLI value for skipping bias correction from `none` to `skip`. Available methods and defaults remain modality-specific.
- DWI and fMRI now use a voxelwise median across 4D volumes to produce the 3D reference, replacing the previous minimum projection.
- Added `--skip-smoothing` to bypass the FSL spatial median filter while still generating the 3D reference image.
- Updated reference and downstream output names to use `MP`, for example `*MP.nii.gz` instead of previous `*DN.nii.gz` references.
- Progress displays use the terminal stream even when logging is enabled.
- FSL BET parameters are printed only when that BET implementation is selected.

## DWI and DSI Studio

- `--denoiser` accepts `patch2self` and `mrtrix3`; denoising remains disabled by default. The batch option is `--dwi-denoiser`.
- MRtrix3 denoises the original 4D DWI before averaging or interpolation, producing `*MR3Denoised.nii.gz` and the noise map `*MR3Noise.nii.gz`. Missing executables, invalid dimensions and external process failures are reported.
- Patch2Self is imported on demand and writes `*P2SDenoised.nii.gz` instead of `*Patch2SelfDenoised.nii.gz`. Headers, float32 data type, qform/sform and units receive explicit handling.
- Added `*P2SNoise.nii.gz`, the signed 4D residual calculated as original minus denoised data. This differs in meaning from the MRtrix noise map.
- b0 averaging resolves gradient files for `P2SDenoised`, `MR3Denoised` and legacy `Patch2SelfDenoised` suffixes, supporting `.bval` and `.btable` inputs.
- DTI registration composes the T2/DTI affine and existing atlas/T2 transformations using `reg_transform -comp`.
- The original Allen annotation and detailed/parental label atlases are resampled directly onto the DTI BET grid using the composed transform. Nearest-neighbor interpolation preserves discrete labels.
- New registration outputs include `*_AtlasToDTI_deformation.nii.gz` and the original unsplit annotation `*_Anno.nii.gz`. The source annotation is checked for existence.
- Improved DSI path resolution. Slice-wise motion correction uses a dedicated temporary directory with more reliable cleanup, replacing global working-directory changes.
- DSI reconstruction passes `--apply_mask`. Isotropic resolution is handled consistently; `auto` uses the same diffusion-derived spacing for reconstruction and ROIs.
- Tracking accepts presets or named `name=value` parameters. Unspecified custom values inherit the `mouse` preset; numeric values, ranges, integer parameters and length limits are validated.
- Current code presets are `default`, `mouse`, `aidamriv2_set`, `test` and `julian_tdcs`. CLI help does not yet list `julian_tdcs`. Removed the former `rat` preset; older preset names and positional parameter lists need review.
- Removed automatic voxel-size-dependent adjustment of tracking step size and minimum length. Resolved tracking values are used directly, and thread count is also passed for the default preset.
- Connectivity now requests `--connectivity_value=all` instead of only `qa` and `count`.
- Stroke lookup tables are created before resampling; label tables accompany resampled seeds. Result files are moved explicitly to the connectivity directory.
- MATLAB-compatible MAT exports use valid, unique variable names, decoded region names and an original-to-exported name mapping. Files are replaced only after successful writing.
- Added `extract_connectivity_metrics.py` to export existing MAT results as cohort CSVs: detailed/parental atlas, `pass`/`end`, region-level `t2r` or region-pair `r2r`, metric selection, optional groups and individual region pairs.
- Region-pair exports exclude diagonals and duplicate directions and validate symmetry. Missing regions remain empty, measured zeros remain zero and missing MAT files are reported.

## Batch Processing

- Added `batchProc.py --exemptionlist` for CSV parameter overrides by `StudyID` and session.
- Validation covers CLI columns, types, boolean values, unknown options, duplicate entries and samples absent from the dataset.
- Unset options use more precise `None` handling to preserve individual script defaults. Help text describes defaults and their sources more explicitly.
- Forwarded new denoising, DWI/fMRI smoothing, fMRI TR and atlas-based T2 masking options.
- The fMRI preprocessing BET selection is no longer automatically forwarded to the subsequent processing stage.
- Explicit CLI options and Git information are logged; available build provenance is copied into the processed project.
- Revised argument groups and CPU/exemption help text.
- Added `batch_DTIdata_extract.py` for FA, MD, RD and AD, including legacy FA0 maps, with metric, ROI, translation table, group and cohort export selection.
- Individual extraction failures do not prevent processing of the remaining cohort. The summary lists failed and skipped entries; failures or entirely unsuccessful runs return an error status.
- Extraction outputs are tracked as a separate stage. Session headers and Git provenance are logged without duplicate child-process entries.

## Incidence and Stroke Mask Analysis

- Added `aggregate_affected_regions.py` to collect existing CSVs from `sub-*/ses-*/anat/affected_Regions/` by filename suffix.
- Exactly one matching file is required per included session. Missing files are reported; ambiguous or invalid data causes errors.
- Output contains subject, session, optional group, total lesion volume, lesion percentage of the brain and consistently sorted regional volume/percentage blocks.
- Regional columns are combined across the cohort; regions missing from individual inputs receive zero. Optional translation tables replace region names with acronyms.
- Groups are assigned per subject from CSV. Duplicate assignments are errors and missing assignments are reported. Existing outputs are replaced only after successful export.
- Added `affected_regions_plots.py` and the aggregation option `--plot` for regional boxplots and lesion-size plots by group/session, exported as PNG and PDF.
- Regional boxplots use logarithmic axes and positive values only; excluded zeros are reported. Lesion-size plots use linear axes and include zeros.
- Boxes show quartiles and medians; whiskers show minima and maxima. Sessions use natural ordering. No significance tests are performed.
- Integrated output tracking and analysis logging into incidence-size extraction. Comments explain the existing incidence registration in more detail.
- Stroke masks can be excluded during DTI ROI analysis and receive special protection during resets.

## T2, T2 Map, fMRI, and ROI Processing

- Added `common/fmri_metadata.py` for shared TR/metadata handling: explicit CLI values take precedence over JSON `RepetitionTime`, with a fallback of 1.42 seconds. TR must be positive and finite.
- fMRI processing supports `--TR`/`--tr`; the batch option is `--func-tr`. Missing metadata or `ScanID` is handled during physiological-data lookup.
- Slice-time correction uses a boolean flag instead of a string parameter; required JSON fields are checked before execution.
- fMRI registration optionally masks the T2 BET using a binarized, slightly dilated atlas annotation: `--atlas-mask-t2`, or `--func-atlas-mask-t2` in batch processing.
- Input matching requires unique results more strictly. External commands are logged and failures propagated. Parental split annotation uses a separate temporary T2 file that is subsequently removed.
- Integrated logging and/or output tracking into the affected T2/T2map registration, extraction, fMRI and ROI scripts.
- `DTIdata_extract.py` writes structured CSVs with `roi_id`, `acronym` and `median`, replacing previous text outputs. Translation-table order and hemisphere-split atlas labels are supported.
- Missing regions and ROIs fully excluded by stroke masks retain empty median cells. Zeros remain measured values; invalid labels and non-finite measurements are checked.
- `DTIdata_extract.py --exclude-stroke-mask PATH` excludes lesion voxels and adds `_exclude-stroke` to the output name. In batch mode, the corresponding flag locates and validates a stroke mask in each DWI folder before extraction.
- Cohort tables combine acronym columns; missing values remain empty. Multiple images of the same metric within a session are treated as ambiguous. Atomic exports respect normal file permissions.
- Extraction logs use `analysis.log` or derived session names; existing processing logs are preserved.
- Timestamps in affected fMRI, ROI, logging and reporting functions use `Europe/Berlin`, including daylight saving time and timezone labels.

## Conversion and Helper Tools

- `conv2Nifti_auto.py` uses `build_conversion_qc_report` from `batch_qc_reports.py`, requesting seven slices instead of ten.
- Removed `plot_sourcedata_niftis.py`. Conversion reports remain available under `Report/Convert2Nifti`.
- Added a post-conversion LIP header check. Filenames containing `fieldmap` or `magnitude` are excluded; orientation mismatches and unreadable files are reported. The check does not reorient images automatically.
- Added parallel processing to `ReorientBatch.py` through `-p/--cpu-percent` (default 50), input validation and orientation summaries.
- Expanded `batch_qc_reports.py` with CLI options, configurable views/report parameters, timestamps, larger images and zoom.
- QC supports corpus callosum views, multiple atlas acronyms, parental/detailed atlas resolution and CSV group filters.
- Added `display_nifti_report.py` for suffix-selected NIfTI reports with optional overlays, transparency, geometry checks and group filters.
- Added `common/artifact_manifest.py` for tracking outputs by modality and processing stage, with `0644` manifest permissions, locking and finalization.
- Added `Reset_proc_folder.py` to plan and apply resets to `base`, `preprocessing`, `registration`, `processing` or `extraction`, including restoration from `brkraw`.
- Resets account for sidecars, gradients and metadata. Untracked files are preserved by default; `--delete-unmanaged` enables explicit cleanup. Stroke masks receive special protection.
- Added `common/script_logging.py` for subject/session log names, terminal output, CLI options, timestamps and Git information.
- Added `common/dti_csv.py` for consistent DTI CSV output naming.

## Breaking or User-visible Changes

- Rebuild the Docker image for the new components and executable paths. Updating the repository does not update an existing container; build provenance identifies the source revision.
- Replace the bias CLI value `none` with `skip`; enable slice-time correction with a flag instead of a string argument.
- Review tracking calls for current presets or named parameters. CLI help does not yet list the implemented `julian_tdcs` preset.
- Median references, revised DTI atlas transformations and tracking parameter handling can change results even with identical raw data.
- Optional denoising, smoothing, TR, masking and stroke-exclusion settings affect results when used.
- Changed `MP`/`P2SDenoised` filenames, additional noise images, CSV formats, log names and expanded connectivity outputs may require downstream script updates.
- Removed `plot_sourcedata_niftis.py`, `iterativeRun_MA.py` and the old example file. Moved three iterative extraction scripts to `legacy/`; the current batch tool is `batch_DTIdata_extract.py`.
- Resets require suitable tracked states or backups and do not replace general data backups.
- Some DSI source comments mention `2026.7.25`, while the Dockerfile downloads the `2025.04.16` CPU build. The Dockerfile determines the built environment.

## Documentation and Figures

- Updated README/manual coverage of Docker provenance, ARM/AMD64 emulation, platform settings and container commands.
- Added documentation for batch logs, denoising, tracking, connectivity CSV extraction and affected-region cohort plots.
- Expanded helper documentation for resets, stroke-mask protection, reorientation, QC CLI options and overlay reports.
- Removed the Gitter reference and updated contact/office-hour information.
- Added update/author information to several script headers. Expanded comments on affine/B-spline registration, transformation direction and resampling.
- No new documentation figures in the net comparison.

## Data and Atlas Files

- `lib/` is unchanged from v3.0; the original release's atlas changes are already part of the baseline.
- Removed previous DTI example data under `bin/3.2.1_DTIdata_extract/testdata/`, including diffusion maps, masks, annotations and text tables.
- Removed tracked macOS `.DS_Store` files from several directories.
- Added `/Publication_AIDAmri_based_on_mouse.pdf` to `.gitignore`; no new PDF in the net comparison.

## Tests

- Added 22 test files since v3.0; Appendix A lists their complete filenames.
- Coverage includes manifests/resets, DWI denoising, DTI registration/extraction, DSI pipeline/batch/MATLAB exports, connectivity, fMRI TR, conversion orientation and QC.
- Additional tests cover affected-region aggregation, plots/statistics, Multiverse, corpus callosum, FA models and upload workflows.
- The presence of tests does not mean every related tool was introduced in this diff; consult the file list for net additions.
- Intermediate ignore-rule changes appear in the history; the final state tracks the test files.
- Tests were not run for this documentation; no runtime validation is claimed.

## Appendix A – Complete File List

Generated with `git diff --name-status --find-renames`. `A` = added, `M` = modified, `D` = deleted, `R100` = moved without content changes.

```text
M	.dockerignore
M	.gitignore
M	ARA/download_ARA.py
M	Dockerfile
M	Manual.md
M	README.md
M	Release_note.md
D	bin/2.1_T2PreProcessing/.DS_Store
M	bin/2.1_T2PreProcessing/MICO.py
M	bin/2.1_T2PreProcessing/applyMICO.py
M	bin/2.1_T2PreProcessing/preProcessing_T2.py
M	bin/2.1_T2PreProcessing/registration_T2.py
M	bin/2.1_T2PreProcessing/t2_value_extraction.py
D	bin/2.2_DTIPreProcessing/.DS_Store
M	bin/2.2_DTIPreProcessing/MICO.py
M	bin/2.2_DTIPreProcessing/applyMICO.py
M	bin/2.2_DTIPreProcessing/averageb0.py
M	bin/2.2_DTIPreProcessing/preProcessing_DTI.py
M	bin/2.2_DTIPreProcessing/registration_DTI.py
D	bin/2.3_fMRIPreProcessing/.DS_Store
M	bin/2.3_fMRIPreProcessing/preProcessing_fMRI.py
M	bin/2.3_fMRIPreProcessing/registration_rsfMRI.py
D	bin/3.1_T2Processing/.DS_Store
A	bin/3.1_T2Processing/affected_regions_plots.py
A	bin/3.1_T2Processing/aggregate_affected_regions.py
M	bin/3.1_T2Processing/getIncidenceSize.py
M	bin/3.1_T2Processing/getIncidenceSize_par.py
D	bin/3.2.1_DTIdata_extract/.DS_Store
M	bin/3.2.1_DTIdata_extract/DTIdata_extract.py
A	bin/3.2.1_DTIdata_extract/batch_DTIdata_extract.py
D	bin/3.2.1_DTIdata_extract/example.txt
D	bin/3.2.1_DTIdata_extract/iterativeRun_MA.py
R100	bin/3.2.1_DTIdata_extract/iterativeRun.py	bin/3.2.1_DTIdata_extract/legacy/iterativeRun.py
R100	bin/3.2.1_DTIdata_extract/iterativeRun_MA_peri-infarct_ROIs.py	bin/3.2.1_DTIdata_extract/legacy/iterativeRun_MA_peri-infarct_ROIs.py
R100	bin/3.2.1_DTIdata_extract/iterativeRun_MA_stroke_mask.py	bin/3.2.1_DTIdata_extract/legacy/iterativeRun_MA_stroke_mask.py
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404.10.1.nii.src.gz.dti.fib.gz.ad.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404.10.1.nii.src.gz.dti.fib.gz.fa0.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404.10.1.nii.src.gz.dti.fib.gz.md.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404.10.1.nii.src.gz.dti.fib.gz.rd.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAllen_scaled.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAnno_rsfMRISplit_scaled.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAnno_rsfMRISplit_scaled.nii.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAnno_rsfMRISplit_scaled_ad.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAnno_rsfMRISplit_scaled_rd.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAnno_rsfMRI_scaled.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAnno_rsfMRI_scaled.nii.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAnno_scaled.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetAnno_scaled.nii.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetMask_scaled.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetStrokeMask_scaled.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetStrokeMask_scaled.nii.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetrsfMRI_Mask_scaled.nii.gz
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetrsfMRI_Mask_scaled.nii.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetrsfMRI_Mask_scaled_ad.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetrsfMRI_Mask_scaled_fa0.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetrsfMRI_Mask_scaled_md.txt
D	bin/3.2.1_DTIdata_extract/testdata/SP_V1_1_1_1_1_20171013_135404SmoothMicoBiasBetrsfMRI_Mask_scaled_rd.txt
D	bin/3.2_DTIConnectivity/.DS_Store
M	bin/3.2_DTIConnectivity/dsi_main.py
M	bin/3.2_DTIConnectivity/dsi_studioPath.txt
M	bin/3.2_DTIConnectivity/dsi_tools.py
A	bin/3.2_DTIConnectivity/extract_connectivity_metrics.py
D	bin/3.3_fMRIActivity/.DS_Store
M	bin/3.3_fMRIActivity/correlate_seed_voxels.py
M	bin/3.3_fMRIActivity/create_seed_rois.py
M	bin/3.3_fMRIActivity/fsl_mean_ts.py
M	bin/3.3_fMRIActivity/getSingleRegTable.py
M	bin/3.3_fMRIActivity/process_fMRI.py
M	bin/3.3_fMRIActivity/regress.py
M	bin/4.1_T2mapPreProcessing/MICO.py
M	bin/4.1_T2mapPreProcessing/applyMICO.py
M	bin/4.1_T2mapPreProcessing/preProcessing_T2MAP.py
M	bin/4.1_T2mapPreProcessing/registration_T2MAP.py
M	bin/4.1_T2mapPreProcessing/t2map_data_extract.py
D	bin/5.1_ROI_analysis/.DS_Store
M	bin/5.1_ROI_analysis/proc_tools.py
D	bin/PV2NIfTiConverter/.DS_Store
D	bin/PV2NIfTiConverter/ReferenceMethods/.DS_Store
M	bin/batchProc.py
A	bin/common/artifact_manifest.py
M	bin/common/bet.py
R100	bin/bet4animal	bin/common/bet4animal
A	bin/common/dti_csv.py
A	bin/common/fmri_metadata.py
A	bin/common/script_logging.py
M	bin/conv2Nifti_auto.py
M	bin/helper_tools/Readme.md
M	bin/helper_tools/ReorientBatch.py
A	bin/helper_tools/Reset_proc_folder.py
M	bin/helper_tools/batch_qc_reports.py
A	bin/helper_tools/display_nifti_report.py
M	bin/helper_tools/getAtlasRegionSize_noBIDS.py
D	bin/helper_tools/plot_sourcedata_niftis.py
M	constraints.txt
A	install/aidamri_entrypoint.sh
M	install/fslinstaller_mod.py
A	install/write_aidamri_git_info.sh
M	requirements.txt
A	tests/test_affected_regions_plots.py
A	tests/test_affected_regions_statistics.py
A	tests/test_aggregate_affected_regions.py
A	tests/test_artifact_manifest.py
A	tests/test_batch_multiverse_qc_reports.py
A	tests/test_batch_qc_reports.py
A	tests/test_cc_fa_group_analysis.py
A	tests/test_conversion_orientation.py
A	tests/test_conversion_qc_reports.py
A	tests/test_create_multiverse_output.py
A	tests/test_display_nifti_report.py
A	tests/test_dsi_batch_integration.py
A	tests/test_dsi_matlab_export.py
A	tests/test_dsi_pipeline.py
A	tests/test_dti_data_extract.py
A	tests/test_dti_registration.py
A	tests/test_dwi_denoising.py
A	tests/test_extract_connectivity_metrics.py
A	tests/test_extract_data_for_multiverse_upload.py
A	tests/test_fa_mixed_models.py
A	tests/test_fmri_tr.py
A	tests/test_iterative_run_ma.py
```

## Appendix B – Complete Commit History

All additional commits, including merges, from oldest to newest. Commit titles may describe intermediate states; the sections above describe the net comparison.

```text
e07c3eb Enhance comments for clarity in incidenceMap function
42e6bc4 Merge branch 'dev-V3'
17f8431 Merge pull request #57 from cramerjulian21c/master
fafa17d feat: add support for custom parameters in QC report generation and update documentation
fbc3e12 docs: update manual and README to include platform warning for ARM users and adjust Docker command
f15036a docs: clarify Docker command for ARM users by specifying platform option
e997499 docs: update README to clarify ARM user warning and container architecture details
2c8c7c4 fix: conditionally print gradient details based on 'bet' argument in preprocessing scripts
3408ed4 Merge pull request #58 from cramerjulian21c/master
8d4944d fix: conditionally print gradient details based on 'bet' argument in preprocessing scripts
8948602 feat: add script logging functionality to preprocessing and registration scripts
5235f0a feat: add timestamp to generated HTML reports in batch QC and NIfTI conversion scripts
9e1ff96 feat: add options to skip smoothing in fMRI and DWI preprocessing
31a691d feat: add option to skip spatial smoothing in DTI preprocessing
b995e69 feat: implement logging functionality in registration scripts
46314de feat: append 'MP' to source_base in DTI and fMRI preprocessing scripts
48471fa feat: implement timezone-aware timestamps in batch QC reports and logging
21b0747 feat: add parallel processing option for batch reorientation of NIfTI files
5503fbb feat: update seaborn version requirement to >=0.13.2 for compatibility
8bed066 feat: add support for exemption lists with parameter overrides in batch processing
b73e76d feat: add Git information tracking and copying in Docker setup
d192857 feat: update seaborn version requirement to 0.13.2 for compatibility
fafe748 Remove tracked DS_Store files
e56b15e feat: update documentation for Docker environment and provenance tracking
5b6574f feat: update date formatting to include month names and set CET timezone for logging
cedf8de feat: update timestamp formatting to use German local time and replace CET with Berlin timezone
4cfb3e8 feat: enhance git info script to extract repository name and remote URL
92c0217 feat: update pillow version to 12.3.0
d142bc6 feat: update exemption list help text to clarify boolean flag usage
fb781a5 feat: update argument handling to allow None checks for configuration options
266cfc1 feat: update slice time correction argument to use boolean flag
3356972 feat: update bias method options to include 'skip' and modify related logic
e9bc8bf feat: add CLI default descriptions and sources for script arguments
0a15820 feat: update isotropic voxel size handling and connectivity values
845bcc5 feat: add apply_mask argument to action=rec command
8549715 feat: update DSI Studio version in Dockerfile and adjust compatibility notes
5918925 feat: update Dockerfile for DSI Studio CPU build and adjust paths
7d316b2 feat: add argument parsing for project-level BET and registration QC reports
5b481b9 feat: update Readme.md to include command-line usage for batch_qc_reports.py
2ee7256 feat: update argument group names in batchProc.py for clarity
84da6c5 Merge pull request #60 from cramerjulian21c/master
d4cc7c7 feat: add artifact manifest tracking and regression tests for AIDAmri processing
88e19fe feat: add Reset_proc_folder.py for manifest-based processing reset and update Readme.md
7f8ed75 feat: integrate output tracking for processing and registration across multiple scripts
8aeaf38 feat: add Publication_AIDAmri_based_on_mouse.pdf to .gitignore
be8baad feat: enhance Reset_proc_folder.py with base restore functionality and update Readme.md
9b81834 feat: update preprocessing scripts to move original NIfTI files and process copies
4c56eab feat: enhance Reset_proc_folder.py to manage metadata files during resets and update Readme.md
5252142 feat: implement --delete-unmanaged option for safe deletion of unmanaged files and directories
5bf3579 feat: enhance Reset_proc_folder.py to preserve stroke mask files during resets
644785b feat: enhance test_artifact_manifest.py to validate preservation of stroke mask files during resets
996b341 feat: update Readme.md to clarify stroke mask preservation and unmanaged file deletion
6c9e40f Merge branch 'Aswendt-Lab:master' into master
192ce14 Merge pull request #61 from cramerjulian21c/master
1cf285e feat: update Readme.md and scripts to clarify preservation of stroke mask files during resets
f74ea3d feat: set manifest file permissions to 644 and add corresponding test
24f153f Merge pull request #63 from cramerjulian21c/master
8e3f2ff feat: log explicit command line options for batch processing
29f5764 feat: log explicit command line options for batch processing
ec0127b feat: dont pass down fmri preprocessing bet setting
05d408b Merge remote-tracking branch 'origin/master'
a20423d feat: enhance logging by adding explicit CLI options and setup for script logging
5f6fbd8 feat: implement dynamic log file naming based on session and subject structure
27fbb36 feat: add support for custom TR input via CLI and enhance metadata handling
b13c929 feat: add option to mask T2 BET with registered atlas annotation during registration
e70a5cd Merge pull request #67 from cramerjulian21c/master
a5145d5 docs: update contact information and clarify office hour details in README
d3711bf docs: remove mention of Gitter Chat from troubleshooting section in README
ce16da3 Merge pull request #70 from cramerjulian21c/master
98f31bc feat: enhance terminal spinner functionality with dynamic stream handling
59e3199 Merge remote-tracking branch 'origin/master'
c9608e6 feat: add corpus callosum QC report generation and update README
61c9afa feat: move bet4animal script to common folder
a008567 feat: add bet4animal script for brain extraction in animal species
e591c44 feat: update regABA2DTI function to include Allen Brain annotation and streamline resampling process
a1f0553 feat: enhance BET_2_MPIreg function with detailed comments for affine and B-spline registration processes
17a4a89 feat: update .gitignore to exclude test directory and remove test_artifact_manifest.py
4a76fba Merge pull request #77 from cramerjulian21c/master
df37d8f feat: enhance DTI preprocessing with support for MRtrix3 dwidenoise and update denoising method options
2a5e311 feat: enhance batch QC reports with image zoom functionality and window-wide images
136644c feat: add git information logging to various scripts for better traceability
4d053ce feat: update DWI denoising method options to support MRtrix3
82d0fb8 feat: enhance Patch2Self denoising to save signed residuals and improve header handling
e12eb82 Merge pull request #80 from cramerjulian21c/master
516a386 feat: implement acronym-based QC report generation and update command-line options
694142f feat: enhance acronym handling in QC report generation to support multiple acronyms
e1a4c33 feat: enhance acronym resolution to support multiple atlas configurations in QC report generation
a9fa35a feat: add display_nifti_report script for generating NIfTI image reports with overlays
f81efde feat: add group filtering functionality to QC report generation and display scripts
606540d Merge pull request #83 from cramerjulian21c/master
bd1b632 feat: add compatibility function for MATLAB in DSI MAT exports
9bc7ccd feat: add script for extracting DSI connectivity metrics to CSV
70c0d55 feat: enhance tracking parameter handling with presets and named values
cfbda21 feat: update tracking parameter help text and remove unused preset rat
67aa8b2 feat: add configuration for julian_tdcs parameters in DSI tools
249e0b6 feat: update documentation to include helper_tools folder for additional scripts and extract_connectivity_metrics explanation
27e9a7d feat: simplify argument parsing in extract_connectivity_metrics script with shorthand options
f03eefc feat: add warnings for missing MAT files and improve CSV export permissions
e1452a6 feat: enhance DTIdata_extract.py with improved error handling and CSV output
2ca12f3 feat: remove obsolete mask and annotation files for rsfMRI analysis
30b1b86 feat: add batch processing script for DTIdata_extract with error handling and CSV output
bf89f70 feat: enhance batch processing and logging in DTIdata_extract scripts with improved error handling and summary reporting
5954015 feat: implement output tracking and error handling in batch extraction scripts
1da5adc feat: enhance batch extraction scripts with session logging and Git information tracking
c7083a5 feat: add cohort CSV output functionality to batch extraction scripts with error handling
6a82128 feat: add group CSV support to batch extraction scripts for cohort output
6e525ab feat: update batch extraction scripts to include extraction stage tracking and logging
9dc263d feat: update argument parser description in batch_DTIdata_extract.py for clarity
27f5c54 feat: remove example.txt and rename iterative run scripts for consistency
4a732e4 feat: update Dockerfile to remove wget command for bet4animal and adjust permissions
de02232 feat: enhance output path handling in batch_DTIdata_extract.py for cohort CSVs
2286f8c feat: update documentation to reference batch_DTIdata_extract.py for iterative processing
804f547 feat: add reference to helper_tools folder in README for additional scripts
2cf7801 feat: add --exclude-stroke-mask option to batch_DTIdata_extract.py and DTIdata_extract.py for stroke exclusion
87d0414 feat: improve handling of missing and fully excluded ROI medians in batch_DTIdata_extract.py and DTIdata_extract.py
68aa07e feat: update logging mechanism to use analysis.log instead of process.log in batch_DTIdata_extract.py and DTIdata_extract.py
d6dc787 feat: add scripts for aggregating and plotting affected regions from CSVs
658aa91 feat: remove tests directory from .gitignore
ae515bc feat: add unit tests for Multiverse registration and CC QC reports, DTI registration, and data extraction for upload
6872c45 Merge pull request #85 from Aswendt-Lab/master
af50778 Merge pull request #86 from Aswendt-Lab/dev-V3
378edb3 Merge branch 'Aswendt-Lab:master' into master
e91832e feat: add latest update information to multiple scripts
c7f5ee2 Merge pull request #87 from cramerjulian21c/master
5fda797 feat: add conversion QC report generation and NIfTI orientation check
4071d1d feat: add unit tests for conversion orientation checks and QC report generation
5d34348 feat: remove plot_sourcedata_niftis.py and update documentation
```
