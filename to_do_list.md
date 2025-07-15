# To Do (Dev)

## Diffusion
 - [ ] Auto-check b-table for best permutation using either DSI Studio's automatic method or by comparing fiber coherence values
   ```python
    # method: 0:DSI, 1:DTI, 4:GQI 7:QSDR, param0: 1.25 (in vivo) diffusion sampling lenth ratio for GQI and QSDR reconstruction, 
    # check_btable: Set –check_btable=1 to test b-table orientation and apply automatic flippin, thread_count: number of multi-threads used to conduct reconstruction
    # flip image orientation in x, y or z direction !! needs to be adjusted according to your data, check fiber tracking result to be anatomically meaningful
    cmd_rec = r'%s --action=%s --source=%s --mask=%s --method=%d --param0=%s --check_btable=%d --half_sphere=%d --cmd=%s'

    # create source files
    filename = os.path.basename(dir_in)
    pos = filename.rfind('.')
    file_src = os.path.join(dir_src, filename[:pos] + ext_src)
    parameters = (dsi_studio, 'src', filename, file_src, b_table)
    os.system(cmd_src % parameters)

    # create fib files
    file_msk = dir_msk
    parameters = (dsi_studio, 'rec', file_src, file_msk, 1, '1.25', 0, 1,'"[Step T2][B-table][flip by]+[Step T2][B-table][flip bz]"')
    os.system(cmd_rec % parameters)
   ```
   Why is the auto-check not preferred? Datasets may vary in which axes to flip in the B-table.
 - [X] Auto populate tracking parameters "min length" to twice the voxel size and "step size" to half voxel size
 - [X] Auto populate `make_isotropic` in `srcgen` for DTI processing
 - [ ] Confirm that tracking options CLI and workflow works as intended
 - [ ] If resampling, the parcellation image does not load for connectivity in DSI Studio. Either resample this or provide a t2 reference image in the connectivity matrix extraction steps.
 - [ ] Image flip Y?
 - [ ] DTI processing does not handle existing outputs well. How should we address this?
 - [ ] FSL slicer images for quality control
 - [ ] DSI Studio fiber images and/or movie for quality control
