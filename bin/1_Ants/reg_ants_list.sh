#!/bin/bash

path_script="$(dirname $0)/reg_ants_20231026.sh"
path_066_anuka="/beegfs/v1/ivnmr_group/studies/066_AnukaMarinaDREADD/processed_data/Resting_state_results"
processed_data="/daten/ivnmr_scratch/michaeld/reg_ants/processed_data"

#$path_script -i ${path_066_anuka}/stroke_cells/2wk/rsfMRI_ANTs $processed_data stroke_cells 2wk AM_30674_1_4_20170820_124922 5 1 8 1
$path_script -i ${path_066_anuka}/stroke_cells/2wk/rsfMRI_ANTs $processed_data stroke_cells 2wk AM_30675_1_3_20170820_135627 5 1 8 1

#$path_script -i ${path_066_anuka}/stroke_cells/6wk/rsfMRI_ANTs $processed_data stroke_cells 6wk AM_30674_1_6_20170918_122017 5 1 8 1
