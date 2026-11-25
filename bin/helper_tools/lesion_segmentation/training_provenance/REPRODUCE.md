# Corrected nnU-Net rerun

This run excludes the user-confirmed no-stroke cases SP_V1_6_2, SP_V1_8_2, SP_V1_10_1. The unresolved SP_T1_4_3 remains withheld. Original validation/test animals and identifiers are preserved. The model starts from fresh random weights, with the same 2D architecture, minibatch, 500-update short-pilot budget, validation checkpoint rule and inference settings as v1. No clicks are used in training or the primary comparison.

From the project directory (preparation refuses to overwrite an existing split):

```sh
.venv/bin/python training/prepare_clean_retrain.py
export MOUSE_STROKE_RUN_DIR="$PWD/training/pilot_v2_clean"
export MOUSE_STROKE_DATASET=Dataset502_MouseStrokeClean
export nnUNet_raw="$MOUSE_STROKE_RUN_DIR/nnunet/raw"
export nnUNet_preprocessed="$MOUSE_STROKE_RUN_DIR/nnunet/preprocessed"
export nnUNet_results="$MOUSE_STROKE_RUN_DIR/nnunet/results"
.venv/bin/nnUNetv2_plan_and_preprocess -d 502 -c 2d --verify_dataset_integrity -np 2 -npfp 2 -gpu_memory_target 4
.venv/bin/python -u training/train_nnunet.py
export MOUSE_STROKE_REVIEW_DIR="$PWD/results/nnunet_refinement_clean_v2"
export MOUSE_STROKE_MODEL_LABEL='nnU-Net 2D clean rerun v2'
.venv/bin/python -u training/prepare_nnunet_refinement.py
.venv/bin/python -u training/compare_clean_nnunet.py
```

The completed v1 dataset and checkpoint are immutable comparison references. New exported training cases, fingerprints, plans, weights and predictions reside under this v2 directory or its separate review output. Do not retrain into a completed run directory.

`comparison/REPORT.md` contains unrefined validation/test comparison first, identical-guidance comparison second. Train-fit and unresolved-case metrics are separate. Containment is intersection/manual reference volume (recall); reverse containment is intersection/prediction volume (precision). Native and mirrored input versions are retained.

This remains a small, single-seed internal pilot with previously reviewed test animals. It is not full-convergence five-fold training or untouched external validation. Sampling/optimization variability and the small cohort limit attribution of score differences to the exclusion alone. Existing GUI predictions are not silently replaced.
