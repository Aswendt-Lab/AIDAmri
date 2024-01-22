#!/bin/bash

# required directory structure and files
# if option -i is specified then files from the input path are copied to the output path
# ${processed_data}/${templates}/
#     $high_res_template
#     $high_res_template_mask
#     $atlas_labels_ambmc
#     $atlas_labels_c57
#     $atlas_to_template_affine
#     $atlas_to_template_warp
# ${processed_data}/${dir_group}/${dir_tp}/${dir_ants}/${dir_seg}/
#     Seg.${high_res_template%_restore.nii*}.nii.gz
#     Seg.Xfm.${subject_anat}.nii.gz
# ${processed_data}/${dir_group}/${dir_tp}/${dir_ants}/${subject}/
#     ${subject_anat}.nii
#     BET.${subject_anat}_restore.nii.gz
#     EPIMean.${subject_func}.nii.gz
#     BET.EPIMean.${subject_func}_restore.nii.gz

dir_tpl='templates'
dir_ants='rsfMRI_ANTs'
dir_seg='_Segmentation'

# template files
high_res_template='MBrTemplate_9T_Cryo_NuNu_8w_n21_256_72_restore.nii'
high_res_template_mask='MBrMask_9T_Cryo_NuNu_8w_n21_256_72.nii.gz'
atlas_labels_ambmc='ambmc-c57bl6-cortex-labels_672_960_rev_lr.nii.gz'
atlas_labels_c57='c57_fixed_labels_resized_672_960_rev.nii.gz'
atlas_to_template_affine='ambmc_rev_to_9T_Cryo_NuNu_8w_0GenericAffine.mat'
atlas_to_template_warp='ambmc_rev_to_9T_Cryo_NuNu_8w_1Warp.nii.gz'

# $0: script filename
script_name="$(basename $0)"

# option parameters and usage
path_in="" # -i
usage0="usage: ./${script_name} [ -i <path_in> ] <processed_data> <dir_group> <dir_tp> <subject> <anat_expno> <anat_procno> <func_expno> <func_procno>"
usage1="-i path_in: path to the input data ANTs directory"
usage2="processed_data: path to the processed data directory"
usage3="dir_group: group subdirectory"
usage4="dir_tp: timepoint subdirectory"
usage5="subject: subject name"
usage6="anat_expno: anatomical data expno"
usage7="anat_procno: anatomical data procno"
usage8="func_expno: functional data expno"
usage9="func_procno: functional data procno"

while getopts ":i:" Option; do
    case $Option in
     i) path_in=$OPTARG;;
    \?) echo "$usage0"
        echo "$usage1"
        echo "$usage2"
        echo "$usage3"
        echo "$usage4"
        echo "$usage5"
        echo "$usage6"
        echo "$usage7"
        echo "$usage8"
        echo "$usage9"
        exit 1;;
    esac
done

# move argument pointer to next index
shift $(($OPTIND - 1))

# check number of arguments
if [ $# -ne 8 ]; then
    echo "$usage0"
    echo "$usage1"
    echo "$usage2"
    echo "$usage3"
    echo "$usage4"
    echo "$usage5"
    echo "$usage6"
    echo "$usage7"
    echo "$usage8"
    echo "$usage9"
    exit 1
fi

echo
if [ "$path_in" != "" ]; then # check if path_in is specified
    echo "./${script_name} -i $path_in $1 $2 $3 $4 $5 $6 $7 $8"
else
    echo "./${script_name} $1 $2 $3 $4 $5 $6 $7 $8"
fi

# $1: processed data directory
processed_data="$1"
if [ ! -d $processed_data ]; then echo "ERROR: directory $processed_data does not exist"; exit 1; fi

# $2: group subdirectory
dir_group="$2"
path_temp="${processed_data}/${dir_group}"
if [ ! -d $path_temp ]; then mkdir $path_temp; fi
if [ ! -d $path_temp ]; then echo "ERROR: directory $path_temp does not exist"; exit 1; fi

# $3: timepoint subdirectory
dir_tp="$3"
path_temp="${processed_data}/${dir_group}/${dir_tp}"
if [ ! -d $path_temp ]; then mkdir $path_temp; fi
if [ ! -d $path_temp ]; then echo "ERROR: directory $path_temp does not exist"; exit 1; fi

path_ants="${processed_data}/${dir_group}/${dir_tp}/${dir_ants}"
if [ ! -d $path_ants ]; then mkdir $path_ants; fi
if [ ! -d $path_ants ]; then echo "ERROR: directory $path_ants does not exist"; exit 1; fi

path_seg="${path_ants}/${dir_seg}"
if [ ! -d $path_seg ]; then mkdir $path_seg; fi
if [ ! -d $path_seg ]; then echo "ERROR: directory $path_seg does not exist"; exit 1; fi

# $4: subject name
subject="$4"
path_sub="${path_ants}/${subject}"
if [ ! -d $path_sub ]; then mkdir $path_sub; fi
if [ ! -d $path_sub ]; then echo "ERROR: directory $path_sub does not exist"; exit 1; fi

# $5: anatomical data expno
anat_expno="$5"

# $6: anatomical data procno
anat_procno="$6"

# $7: functional data expno
func_expno="$7"

# $8: functional data procno
func_procno="$8"

subject_anat="${subject}.${anat_expno}.${anat_procno}"
subject_func="${subject}.${func_expno}.${func_procno}"

# segmented files
seg_fixfile="Seg.${high_res_template%_restore.nii*}.nii.gz"
seg_movfile="Seg.Xfm.${subject_anat}.nii.gz"

# subject files
sub_anat="${subject_anat}.nii"
sub_anat_bet_bc="BET.${subject_anat}_restore.nii.gz"
sub_func_mean="EPIMean.${subject_func}.nii.gz"
sub_func_mean_bet_bc="BET.EPIMean.${subject_func}_restore.nii.gz"

# copy input files to output directory
if [ "$path_in" != "" ]; then # check if path_in is specified
    if [ -f ${path_in}/${dir_seg}/${seg_fixfile} ] && [ ! -f ${path_seg}/${seg_fixfile} ]; then cp ${path_in}/${dir_seg}/${seg_fixfile} $path_seg; fi
    if [ -f ${path_in}/${dir_seg}/${seg_movfile} ] && [ ! -f ${path_seg}/${seg_movfile} ]; then cp ${path_in}/${dir_seg}/${seg_movfile} $path_seg; fi
    if [ -f ${path_in}/${subject}/${sub_anat} ] && [ ! -f ${path_sub}/${sub_anat} ]; then cp ${path_in}/${subject}/${sub_anat} $path_sub; fi
    if [ -f ${path_in}/${subject}/${sub_anat_bet_bc} ] && [ ! -f ${path_sub}/${sub_anat_bet_bc} ]; then cp ${path_in}/${subject}/${sub_anat_bet_bc} $path_sub; fi
    if [ -f ${path_in}/${subject}/${sub_func_mean} ] && [ ! -f ${path_sub}/${sub_func_mean} ]; then cp ${path_in}/${subject}/${sub_func_mean} $path_sub; fi
    if [ -f ${path_in}/${subject}/${sub_func_mean_bet_bc} ] && [ ! -f ${path_sub}/${sub_func_mean_bet_bc} ]; then cp ${path_in}/${subject}/${sub_func_mean_bet_bc} $path_sub; fi
    sleep 5 # wait 5 seconds
fi

path_tpl="${processed_data}/${dir_tpl}"

# check if template files do exist
if [ ! -f ${path_tpl}/${high_res_template} ]; then echo "ERROR: template file ${path_tpl}/${high_res_template} does not exist"; exit 1; fi
if [ ! -f ${path_tpl}/${high_res_template_mask} ]; then echo "ERROR: template file ${path_tpl}/${high_res_template_mask} does not exist"; exit 1; fi
if [ ! -f ${path_tpl}/${atlas_labels_ambmc} ]; then echo "ERROR: template file ${path_tpl}/${atlas_labels_ambmc} does not exist"; exit 1; fi
if [ ! -f ${path_tpl}/${atlas_labels_c57} ]; then echo "ERROR: template file ${path_tpl}/${atlas_labels_c57} does not exist"; exit 1; fi
if [ ! -f ${path_tpl}/${atlas_to_template_affine} ]; then echo "ERROR: template file ${path_tpl}/${atlas_to_template_affine} does not exist"; exit 1; fi
if [ ! -f ${path_tpl}/${atlas_to_template_warp} ]; then echo "ERROR: template file ${path_tpl}/${atlas_to_template_warp} does not exist"; exit 1; fi

# check if segmented files do exist
if [ ! -f ${path_seg}/${seg_fixfile} ]; then echo "ERROR: segmented file ${path_seg}/${seg_fixfile} does not exist"; exit 1; fi
if [ ! -f ${path_seg}/${seg_movfile} ]; then echo "ERROR: segmented file ${path_seg}/${seg_movfile} does not exist"; exit 1; fi

# check if subject files do exist
if [ ! -f ${path_sub}/${sub_anat} ]; then echo "ERROR: subject file ${path_sub}/${sub_anat} does not exist"; exit 1; fi
if [ ! -f ${path_sub}/${sub_anat_bet_bc} ]; then echo "ERROR: subject file ${path_sub}/${sub_anat_bet_bc} does not exist"; exit 1; fi
if [ ! -f ${path_sub}/${sub_func_mean} ]; then echo "ERROR: subject file ${path_sub}/${sub_func_mean} does not exist"; exit 1; fi
if [ ! -f ${path_sub}/${sub_func_mean_bet_bc} ]; then echo "ERROR: subject file ${path_sub}/${sub_func_mean_bet_bc} does not exist"; exit 1; fi

outbase_rigid="${dir_group}_${dir_tp}_ants_3d_rigid_"
outbase_affine="${dir_group}_${dir_tp}_ants_3d_affine_"
outbase_pse="${dir_group}_${dir_tp}_ants_3d_pse_"
outbase_diff="${dir_group}_${dir_tp}_ants_3d_diff_"

cd $path_sub

# 1. rsfMRI EPI mean -> anatomical T2w (rigid)
antsRegistrationSyN.sh -n 4 -t r -d 3 -p f -f $sub_anat_bet_bc -m $sub_func_mean_bet_bc -o $outbase_rigid > ${outbase_rigid}log.txt
# apply transformation
${ANTSPATH}/antsApplyTransforms -d 3 -e 0 --float 1 -i $sub_func_mean -o ${outbase_rigid}Xfm.nii.gz -r $sub_anat -t ${outbase_rigid}0GenericAffine.mat

# 2. a) anatomical T2w -> in-house created MBrTemplate (rigid + affine)
antsRegistrationSyN.sh -n 4 -t a -d 3 -f ${path_tpl}/${high_res_template} -m $sub_anat_bet_bc -o $outbase_affine > ${outbase_affine}log.txt
# apply transformation
${ANTSPATH}/antsApplyTransforms -d 3 -e 0 --float 1 -i $sub_anat -o ${outbase_affine}Xfm.nii.gz -r ${path_tpl}/${high_res_template} -n Linear -t ${outbase_affine}0GenericAffine.mat

# overwrite header of ${outbase_affine}Xfm
d1=$(fslval ${outbase_affine}Xfm dim1); d2=$(fslval ${outbase_affine}Xfm dim2); d3=$(fslval ${outbase_affine}Xfm dim3); dt=$(fslval ${outbase_affine}Xfm datatype)
pd1=$(fslval ${outbase_affine}Xfm pixdim1); pd2=$(fslval ${outbase_affine}Xfm pixdim2); pd3=$(fslval ${outbase_affine}Xfm pixdim3); pd4=$(fslval ${outbase_affine}Xfm pixdim4)
fslcreatehd $d1 $d2 $d3 1 $pd1 $pd2 $pd3 $pd4 0 0 0 $dt ${outbase_affine}Xfm
# apply template mask
fslmaths ${outbase_affine}Xfm -mas ${path_tpl}/${high_res_template_mask} ${outbase_affine}Xfm_Msk
# bias correction output ${outbase_affine}restore
fast -l 20 -t 2 --nopve -B -o ${outbase_affine%_} ${outbase_affine}Xfm_Msk
rm ${outbase_affine}seg.nii.gz

#fixfile="../${dir_seg}/${seg_fixfile}"
fixfile="${path_seg}/${seg_fixfile}"
echo "fixfile: $fixfile"

#movfile="../${dir_seg}/${seg_movfile}"
movfile="${path_seg}/${seg_movfile}"
echo "movfile: $movfile"

# 2. b) segmented transformed anatomical T2w -> segmented in-house created MBrTemplate (PSE)
${ANTSPATH}/ANTS 3 -o $outbase_pse -i 91x70x55x40x30 -r Gauss[3,0.5] -t SyN[0.2] -m PSE[${fixfile},${movfile},${fixfile},${movfile},0.75,0.1,11,0,10] -m MSQ[${fixfile},${movfile},3,0] --number-of-affine-iterations 0 > ${outbase_pse}log.txt
# transformations
${ANTSPATH}/WarpImageMultiTransform 3 $movfile ${outbase_pse}Warped.nii.gz -R $fixfile ${outbase_pse}Warp.nii.gz ${outbase_pse}Affine.txt >> ${outbase_pse}log.txt
${ANTSPATH}/WarpImageMultiTransform 3 $fixfile ${outbase_pse}InverseWarped.nii.gz -R $movfile -i ${outbase_pse}Affine.txt ${outbase_pse}InverseWarp.nii.gz >> ${outbase_pse}log.txt
${ANTSPATH}/WarpImageMultiTransform 3 ${outbase_affine}restore.nii.gz ${outbase_pse}Xfm_Warped.nii.gz -R ${path_tpl}/${high_res_template} ${outbase_pse}Warp.nii.gz ${outbase_pse}Affine.txt >> ${outbase_pse}log.txt
${ANTSPATH}/WarpImageMultiTransform 3 ${path_tpl}/${high_res_template} ${outbase_pse}Xfm_InverseWarped.nii.gz -R ${outbase_affine}restore.nii.gz -i ${outbase_pse}Affine.txt ${outbase_pse}InverseWarp.nii.gz >> ${outbase_pse}log.txt
# create grid for visualization
${ANTSPATH}/CreateWarpedGridImage 3 ${outbase_pse}Warp.nii.gz ${outbase_pse}grid.nii.gz
${ANTSPATH}/CreateWarpedGridImage 3 ${outbase_pse}InverseWarp.nii.gz ${outbase_pse}grid_inv.nii.gz

movfile=${outbase_pse}Xfm_Warped.nii.gz
echo "movfile: $movfile"

# 2. c) transformed anatomical T2w -> in-house created MBrTemplate (diffeomorphic)
antsRegistrationSyN.sh -n 4 -t so -d 3 -f ${path_tpl}/${high_res_template} -m $movfile -o $outbase_diff > ${outbase_diff}log.txt
# create grid for visualization
${ANTSPATH}/CreateWarpedGridImage 3 ${outbase_diff}1Warp.nii.gz ${outbase_diff}grid.nii.gz
${ANTSPATH}/CreateWarpedGridImage 3 ${outbase_diff}1InverseWarp.nii.gz ${outbase_diff}grid_inv.nii.gz

# apply transformations to anatomical T2w -> in-house created MBrTemplate
${ANTSPATH}/antsApplyTransforms -d 3 -e 0 --float 1 -i $sub_anat -o Warped_${sub_anat} -r ${path_tpl}/${high_res_template} -n Linear -t ${outbase_diff}1Warp.nii.gz -t ${outbase_diff}0GenericAffine.mat -t ${outbase_pse}Warp.nii.gz -t ${outbase_affine}0GenericAffine.mat

# apply transformations to atlas labels -> anatomical T2w
${ANTSPATH}/antsApplyTransforms -d 3 -e 0 --float 1 -i ${path_tpl}/${atlas_labels_ambmc} -o Warped_${atlas_labels_ambmc} -r $sub_anat -n NearestNeighbor -t [${outbase_affine}0GenericAffine.mat,1] -t ${outbase_pse}InverseWarp.nii.gz -t [${outbase_diff}0GenericAffine.mat,1] -t ${outbase_diff}1InverseWarp.nii.gz -t ${path_tpl}/${atlas_to_template_warp} -t ${path_tpl}/${atlas_to_template_affine}
${ANTSPATH}/antsApplyTransforms -d 3 -e 0 --float 1 -i ${path_tpl}/${atlas_labels_c57} -o Warped_${atlas_labels_c57} -r $sub_anat -n NearestNeighbor -t [${outbase_affine}0GenericAffine.mat,1] -t ${outbase_pse}InverseWarp.nii.gz -t [${outbase_diff}0GenericAffine.mat,1] -t ${outbase_diff}1InverseWarp.nii.gz -t ${path_tpl}/${atlas_to_template_warp} -t ${path_tpl}/${atlas_to_template_affine}

# apply transformations to atlas labels -> rsfMRI EPI
${ANTSPATH}/antsApplyTransforms -d 3 -e 0 --float 1 -i ${path_tpl}/${atlas_labels_ambmc} -o Warped.cortex.reg.${subject_func}.nii.gz -r $sub_func_mean -n NearestNeighbor -t [${outbase_rigid}0GenericAffine.mat,1] -t [${outbase_affine}0GenericAffine.mat,1] -t ${outbase_pse}InverseWarp.nii.gz -t [${outbase_diff}0GenericAffine.mat,1] -t ${outbase_diff}1InverseWarp.nii.gz -t ${path_tpl}/${atlas_to_template_warp} -t ${path_tpl}/${atlas_to_template_affine}
${ANTSPATH}/antsApplyTransforms -d 3 -e 0 --float 1 -i ${path_tpl}/${atlas_labels_c57} -o Warped.striatum.reg.${subject_func}.nii.gz -r $sub_func_mean -n NearestNeighbor -t [${outbase_rigid}0GenericAffine.mat,1] -t [${outbase_affine}0GenericAffine.mat,1] -t ${outbase_pse}InverseWarp.nii.gz -t [${outbase_diff}0GenericAffine.mat,1] -t ${outbase_diff}1InverseWarp.nii.gz -t ${path_tpl}/${atlas_to_template_warp} -t ${path_tpl}/${atlas_to_template_affine}