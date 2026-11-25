# Mouse lesion guidance — version 1.2.0 research pilot

Part of [AIDAmri](../../../README.md), under `bin/helper_tools/lesion_segmentation`.
After cloning or downloading AIDAmri, open this folder in Finder. For the guided
demo, double-click **Setup-demo.command**, then **Start.command**. For new-image
detection, use **Setup.command** to install the full inference environment.
All example data, saved guidance and the frozen checkpoint are included here;
no separate archive or model download is required.

This local browser application lets you place stroke seeds and avoidance ROIs,
then run a frozen clean nnU-Net and select a continuous lesion region. Everything
stays on your computer. No account or cloud service is needed.

## Start with the saved examples

1. Extract the whole folder to a writable location. Keep its folders together.
2. Install Python 3.10, 3.11 or 3.12 if needed. On macOS, double-click
   **Setup.command** once. It creates an environment and installs the
   review and detection dependencies (internet needed for this setup).
3. Double-click **Start.command**. Your browser opens http://127.0.0.1:8780.
   Leave the package folder in place while the application is running.
4. Click **Start guided demo** to learn with an included example, or click
   **View example result** to inspect the original saved result immediately.
   For normal work, choose a dataset. Place a cyan lesion seed inside the stroke. Select
   Avoidance ROI to mark unwanted areas, such as ventricles. Its radius is in
   millimetres. Navigate slices and add more marks as needed.
5. Click **Run guided segmentation + save**. Both input directions are saved.
   Check seed-hit counts, Dice and containment. Manual-reference and original
   baseline contours can be toggled on; the reference starts hidden.
6. To replay the supplied human guidance, choose a saved run, click Restore
   guidance for editing, then Run. Editing marks hides the earlier result until
   another explicit run. New runs never overwrite the supplied sessions.

Terminal alternative (macOS/Linux; use the equivalent Python executable on Windows):

```sh
python3.11 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/start_validation_guidance.py --dataset examples --port 8780
```

The examples already contain frozen network predictions. Reviewing them does not
require PyTorch or a GPU. The launcher runs the local server detached; closing the
browser leaves it available. Its PID and logs are in `.gui_runtime/`. To stop it,
use your process manager or `kill` with the PID from `validation_8780.pid`.
If another dataset/server uses this port, start with a different `--port` value.

## Guided demo — learn without loading your own files

For a small demo-only installation, double-click **Setup-demo.command** instead
of Setup.command. This installs only NumPy, SciPy and NiBabel; neither PyTorch nor
a GPU is needed for the included example. Then open Start.command.

Click **Start guided demo** and follow the five steps:

1. Inspect the included T2w MRI with the slice controls.
2. Place a lesion seed inside the dashed cyan circle, or choose **Use recorded
   seed** to copy the original user's mark.
3. Add an avoidance sphere at the dashed orange target, or choose **Use recorded
   avoidance ROI**. The physical radius is shown in millimetres.
4. Click **Run guided segmentation + save**. Both input directions are processed
   from the example's cached frozen predictions. Practice output goes into
   `examples/guided_continuity/demo_sessions/`, separate from the research review.
5. Compare reference/baseline/guided contours and Dice/containment. The recorded
   human marks yield native Dice **0.380 → 0.851** and containment **89.9%** on
   GV_T3_3_2_1_2_20181121_091427. Different practice clicks may give other scores.

**View example result** / **View recorded result** loads the actual preserved
session immediately, with both input directions and all four metric rows. This
is read only: no new run or model inference occurs. **Try yourself / restart**
clears tutorial marks only. **Exit demo** restores the prior scan, its draft marks,
slice and result. Replaying or practising never overwrites recorded guidance.
The example is previously reviewed learning data, not independent validation.

All 11 reviewed scans, their manual masks, baseline predictions and 12 original
human sessions remain included. The guided tour uses one reproducible example;
other examples can be explored through the dataset and saved-session selectors.
Run Setup.command later to install full detection dependencies for new scans.

## New scans — no command line needed

1. Click **Load new scan…** in the interface.
2. Choose a skull-stripped T2-weighted NIfTI image (`.nii` or `.nii.gz`).
3. Optionally choose the matching manual stroke mask; leave it empty to segment
   without a reference. Click **Load image**.
4. The new scan appears in the dataset selector. Place positive lesion seeds
   and avoidance ROIs, then click **Run guided segmentation + save**.

Image and optional mask are copied locally into `examples/imported_scans/`.
New scans persist after reopening the app and never replace frozen examples.
Geometry/data errors are shown in the dialog; a rejected upload adds no scan.
Files may be up to 128 MB each and 32 million voxels. A supplied reference must
already match MRI shape and affine. No automatic mask resampling is performed.
With no manual reference, overlap scores display N/A.

No detection runs while loading files or placing seeds. After Run, the checkpoint
is applied if predictions do not yet exist, then guidance selects the region.
First inference can take longer on CPU. Apple MPS, CUDA and CPU are selected
in that order of availability (CUDA preferred when available). macOS/Apple GPU
execution has been tested; other devices may have small boundary differences.

Setup.command installs inference dependencies as well as review dependencies.
If you used the lightweight terminal installation above or are upgrading an
older package, run Setup.command to install the full environment before loading
new scans. A missing dependency is explained in the interface and guidance is
saved before detection is attempted.

Use the same skull-stripped T2 preparation and anatomical slice-stack convention
as the examples (native voxel axis 2). Raw scanner/DTI files and arbitrary
reorientation are not supported. Direction harmonization remains deferred.
Both native and axis-0-mirrored versions are retained; do not select a test
orientation using its reference score.

The optional `scripts/import_guidance_dataset.py` command remains available for
batch imports, but is not needed for the file picker workflow.

## What the guidance does

Avoidance spheres are cut out first. The first positive seed that intersects the
network prediction anchors a continuous 3D track. Additional reachable seeds guide
branch selection. There is at most one 8-connected region per native slice, and
neighbouring regions must touch under 26-connectivity. No jumping across an empty
slice or to a distant location, and no snapping of missed seeds, is allowed.
Avoidance-only runs are possible but unanchored; the largest continuous track
can be the wrong region, so a positive seed is recommended.

This is manual selection from a frozen prediction, not seed-conditioned neural
inference or retraining. It cannot add tissue the model missed. It can stop early
at a prediction gap. Incompatible/missed seeds produce explicit warnings.

## Outputs and saved status

Every Run writes a new `guided_continuity/manual_sessions/<time_id>/` inside the
selected dataset: guidance JSON, native/mirrored NIfTI masks, avoidance mask,
baseline/guided metric CSV, checkpoint/algorithm provenance and warnings.
Guidance is saved before any new inference. Inputs and reference masks remain
unchanged. Browser draft marks persist locally; Run saves them to disk.

The frozen `review/` report uses the latest of 12 human sessions for each of
11 P1 scans (one retry). Mean native Dice: **0.454 baseline → 0.836 guided**;
guided containment: **93.2%**; prediction containment: **77.4%**.
Eight reviewed scans belong to the earlier strict-new subset: Dice **0.462 →
0.823**, containment **93.6%**. The other three overlap development data.
All 24 saved masks and all 48 metric rows were verified. One early failed seed
and its retry remain preserved. These are user-selected, manually refined P1
examples, not independent validation; P7 manual review remains unperformed.

Containment = intersection/manual lesion volume (recall). Reverse containment =
intersection/predicted volume (precision). Undefined denominators show N/A.
The user's assessment “working perfectly” is preserved in STATUS.json alongside
the measured scores and their scope. This remains a short-trained research pilot.

## Reproduce and verify

```sh
.venv/bin/python scripts/reproduce_guidance_examples.py
```

This checks the release checksums and reproduces all supplied masks from saved
guidance, including the failed first attempt, and rechecks all metrics. It does
not create new human sessions or modify the examples. To also rerun the network
on all 11 example MRI volumes (after installing inference dependencies):

```sh
.venv/bin/python scripts/reproduce_guidance_examples.py --network
```

The network rerun uses temporary output and writes `network_reproduction.json`;
it never overwrites the baseline. Cross-device floating-point differences may
affect boundary voxels. The checkpoint, nnU-Net plans/dataset metadata, source
images, manual masks, baseline predictions, all human sessions, code, dependency
versions and training provenance are bundled. Full retraining is not included;
this handoff reproduces frozen-model inference and manual refinement.

See `PROVENANCE.json`, `training_provenance/`, `SHA256SUMS.json` and
`THIRD_PARTY_NOTICES.md`. Do not remove the example inputs needed for reproduction.

## Publication provenance

Historical training/audit metadata use portable `historical_project/` and
`historical_training_inputs/` labels for original development material that is
not bundled. `review/source_checksums.json` records the original pre-publication
source audit; JSON metadata paths were made portable for this repository.
`SHA256SUMS.json` verifies the actual published files. The MRI, reference,
prediction and model arrays are unchanged. This directory is independent of the
rest of the AIDAmri processing environment. Generated uploads, demo practice runs,
new sessions, environments and runtime caches are ignored by Git.
