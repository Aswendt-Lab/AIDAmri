#!/usr/bin/env python3
"""Run per-slice motion and atlas-based CompCor regression on processed rat EPI."""

import argparse
import json
import re
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.linalg import svd
from scipy.signal import butter, detrend, sosfiltfilt


WHITE_MATTER_LABELS = (781, 782, 891, 892, 941, 942, 1101, 1102)
VENTRICULAR_LABELS = (1171, 1172)
SEED_LABELS = {"S1_L": 521, "S1_R": 522, "ACA_R": 502}
MOTION_COLUMNS = (2, 3, 4)  # MCFLIRT Rz, Tx, Ty
PASSBAND_HZ = (0.01, 0.20)
DISCARDED_VOLUMES = 5


def data_files(directory, pattern):
    return sorted(path for path in directory.glob(pattern)
                  if not path.name.startswith("._") and path.is_file())


def one_file(directory, pattern):
    matches = data_files(directory, pattern)
    if len(matches) != 1:
        raise ValueError(f"Expected one {pattern} in {directory}, found {len(matches)}")
    return matches[0]


def spatial_data(path, reference):
    image = nib.load(path)
    if image.shape != reference.shape[:3]:
        raise ValueError(f"Spatial shape mismatch: {path}")
    corners = np.array(np.meshgrid(*[(0, n - 1) for n in image.shape])).reshape(3, -1).T
    displacement = np.linalg.norm(
        nib.affines.apply_affine(image.affine, corners)
        - nib.affines.apply_affine(reference.affine, corners), axis=1
    )
    if displacement.max() > 0.1 * min(reference.header.get_zooms()[:3]):
        raise ValueError(f"Spatial affine mismatch: {path}")
    return np.asanyarray(image.dataobj)


def get_inputs(func):
    native = func / "rs-fMRI_niiData"
    # Prefer the slice-time-corrected, trimmed native input used by the exporter.
    pattern = "*_EPI_mcf_st_f.nii.gz"
    if not data_files(native, pattern):
        pattern = "*_EPI_mcf_f.nii.gz"
    path = one_file(native, pattern)
    stem = path.name.removesuffix(".nii.gz")
    raw_stem = stem.split("_mcf", 1)[0]
    # Reuse the exporter's masked copy when it exists, reducing volume reads.
    masked = native / "Multiverse_Output" / f"{stem}_BET.nii.gz"
    return (
        masked if masked.is_file() else path,
        one_file(func, f"{raw_stem}*Bet_mask.nii.gz"),
        one_file(func, f"{raw_stem}*Bet_Anno.nii.gz"),
        func / "rs-fMRI_mcf" / f"{raw_stem}_mcf.mat",
        func / f"{raw_stem}.json",
    )


def repetition_time(image, sidecar):
    header_tr = float(image.header.get_zooms()[3])
    unit = image.header.get_xyzt_units()[1]
    if unit == "sec":
        return header_tr
    if unit == "msec":
        return header_tr / 1000

    # Some intermediate FSL outputs lose temporal units. Use the BIDS sidecar
    # as the authoritative source instead of assuming the header is in seconds.
    metadata = json.loads(sidecar.read_text())
    if "RepetitionTime" in metadata:
        return float(metadata["RepetitionTime"])
    timing = np.asarray(metadata.get("VolumeTiming", []), dtype=float)
    if timing.shape != (image.shape[3] + DISCARDED_VOLUMES,):
        raise ValueError(f"No usable repetition time in {sidecar}")
    intervals = np.diff(timing)
    if not np.allclose(intervals, intervals[0], atol=1e-4, rtol=1e-4):
        raise ValueError(f"Nonuniform volume timing in {sidecar}")
    return float(intervals[0])


def read_motion(directory, n_slices, n_volumes):
    by_slice = {}
    for path in data_files(directory, "*.par"):
        match = re.search(r"(\d{4})\.nii(?:\.gz)?\.par$", path.name)
        if match is None or int(match[1]) in by_slice:
            raise ValueError(f"Unexpected or duplicate slice motion filename: {path}")
        by_slice[int(match[1])] = path
    if set(by_slice) != set(range(n_slices)):
        raise ValueError(f"Motion files do not cover all {n_slices} native slices: {directory}")
    traces = []
    for index in range(n_slices):
        table = np.loadtxt(by_slice[index], ndmin=2)
        if table.shape != (n_volumes + DISCARDED_VOLUMES, 6) or not np.isfinite(table).all():
            raise ValueError(f"Unexpected motion table shape or values: {by_slice[index]}")
        traces.append(table[DISCARDED_VOLUMES:, MOTION_COLUMNS])
    return traces


def components(series, count):
    sd = series.std(axis=0)
    columns = sd > np.finfo(np.float32).eps
    if columns.sum() < count:
        raise ValueError("Too few variable WM/CSF time courses")
    standardized = series[:, columns] / sd[columns]
    left, singular, _ = svd(standardized, full_matrices=False)
    if singular[count - 1] <= singular[0] * max(standardized.shape) * np.finfo(singular.dtype).eps:
        raise ValueError("WM/CSF time courses have insufficient rank")
    return left[:, :count]


def subtract_nuisance(data, regressors, ridge):
    centered = regressors - regressors.mean(axis=0)
    sd = centered.std(axis=0)
    valid = sd > 1e-10
    design = centered[:, valid] / sd[valid]
    if ridge:
        penalty = ridge * np.eye(design.shape[1])
        weights = np.linalg.solve(design.T @ design + penalty, design.T @ data)
    else:
        weights = np.linalg.lstsq(design, data, rcond=1e-8)[0]
    return data - design @ weights


def pair_correlations(series):
    results = {}
    for a, b in (("S1_L", "S1_R"), ("S1_R", "ACA_R")):
        if a in series and b in series and min(series[a].std(), series[b].std()) > 0:
            results[f"{a}__{b}"] = float(np.corrcoef(series[a], series[b])[0, 1])
    return results


def run_scan(func, output, n_components, ridge):
    source, mask_file, atlas_file, motion_dir, sidecar = get_inputs(func)
    image = nib.load(source)
    if len(image.shape) != 4:
        raise ValueError(f"Expected a 4D native EPI: {source}")
    tr = repetition_time(image, sidecar)
    if not np.isfinite(tr) or tr <= 0 or PASSBAND_HZ[1] >= 0.5 / tr:
        raise ValueError(f"Invalid TR or passband exceeds Nyquist at TR={tr:g} s: {source}")
    edge = round(30 / tr)
    if image.shape[3] <= 2 * edge + 5:
        raise ValueError(f"Too few frames after trimming filter edges: {source}")
    filter_sos = butter(2, PASSBAND_HZ, btype="bandpass", fs=1 / tr, output="sos")

    def bandpass(values):
        result = sosfiltfilt(filter_sos, detrend(values, axis=0), axis=0)
        return result - result.mean(axis=0)

    mask = spatial_data(mask_file, image) > 0
    labels = spatial_data(atlas_file, image)
    data = image.get_fdata(dtype=np.float32, caching="unchanged")
    if not np.isfinite(data[mask]).all():
        raise ValueError(f"Non-finite values in the brain mask: {source}")
    mask &= (data.max(axis=-1) - data.min(axis=-1)) > 0
    wm = mask & np.isin(labels, WHITE_MATTER_LABELS)
    csf = mask & np.isin(labels, VENTRICULAR_LABELS)
    if not np.any(wm) or not np.any(csf):
        raise ValueError(f"Atlas and EPI coverage leave an empty WM or CSF mask: {source}")
    compcor = components(bandpass(data[wm | csf].T), n_components)
    motion = [
        bandpass(trace) for trace in read_motion(motion_dir, image.shape[2], image.shape[3])
    ]

    corrected = np.zeros_like(data)
    before_seeds, after_seeds = {}, {}
    before_sum = np.zeros((len(SEED_LABELS), image.shape[3]))
    after_sum = np.zeros_like(before_sum)
    seed_counts = np.zeros(len(SEED_LABELS), dtype=int)
    removed_variance = []
    for index, trace in enumerate(motion):
        active = mask[:, :, index]
        if not np.any(active):
            continue
        original = data[:, :, index][active].T
        filtered = bandpass(original)
        residual = subtract_nuisance(filtered, np.column_stack([trace, compcor]), ridge)
        corrected[:, :, index][active] = (residual + original.mean(axis=0)).T

        power = np.square(filtered).sum(axis=0)
        variable = power > 0
        removed_variance.extend(
            1 - np.square(residual[:, variable]).sum(axis=0) / power[variable]
        )
        slice_labels = labels[:, :, index][active]
        for seed_index, label in enumerate(SEED_LABELS.values()):
            selected = slice_labels == label
            seed_counts[seed_index] += selected.sum()
            before_sum[seed_index] += filtered[:, selected].sum(axis=1)
            after_sum[seed_index] += residual[:, selected].sum(axis=1)

    # Omit filter edges from the descriptive connectivity metrics.
    for index, name in enumerate(SEED_LABELS):
        if seed_counts[index]:
            before_seeds[name] = before_sum[index, edge:-edge] / seed_counts[index]
            after_seeds[name] = after_sum[index, edge:-edge] / seed_counts[index]

    qc = {
        "source": str(source),
        "brain_mask": str(mask_file),
        "tissue_atlas": str(atlas_file),
        "motion_directory": str(motion_dir),
        "tr_seconds": tr,
        "volumes": image.shape[3],
        "passband_hz": PASSBAND_HZ,
        "aCompCor_components": n_components,
        "ridge_penalty": ridge,
        "brain_voxels": int(mask.sum()),
        "white_matter_voxels": int(wm.sum()),
        "csf_voxels": int(csf.sum()),
        "seed_voxels": {name: int(n) for name, n in zip(SEED_LABELS, seed_counts)},
        "median_in_sample_voxel_variance_removed": float(np.median(removed_variance)),
        "bandpass_only_correlations": pair_correlations(before_seeds),
        "denoised_correlations": pair_correlations(after_seeds),
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".regression-", dir=output.parent) as staging:
        stage = Path(staging)
        denoised = nib.Nifti1Image(corrected, image.affine, header=image.header.copy())
        denoised.set_data_dtype(np.float32)
        denoised.header.set_zooms((*image.header.get_zooms()[:3], tr))
        denoised.header.set_xyzt_units(image.header.get_xyzt_units()[0], "sec")
        denoised.set_sform(*image.get_sform(coded=True))
        denoised.set_qform(*image.get_qform(coded=True))
        stem = source.name.removesuffix(".nii.gz").removesuffix("_BET")
        nib.save(denoised, stage / f"{stem}_nuisance-regressed_bold.nii.gz")
        columns = [f"aCompCor_{i + 1}" for i in range(n_components)]
        columns += [f"slice_{i:02d}_{p}" for i in range(len(motion)) for p in ("Rz", "Tx", "Ty")]
        np.savetxt(stage / "nuisance_regressors.tsv", np.column_stack([compcor, *motion]),
                   delimiter="\t", fmt="%.9g", header="\t".join(columns), comments="")
        (stage / "regression_qc.json").write_text(json.dumps(qc, indent=2, allow_nan=False) + "\n")
        # Publish a complete subject/session result atomically.
        stage.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path, help="Processed dataset containing sub-*/ses-*/func")
    parser.add_argument("output_root", type=Path, help="Separate directory for corrected results")
    parser.add_argument("--subject", action="append", help="Process a subject ID, e.g. sub-300100; repeatable")
    parser.add_argument("--components", type=int, default=3, help="Number of WM/CSF PCs (default: 3)")
    parser.add_argument("--ridge", type=float, default=0.0,
                        help="Penalty on standardized nuisance coefficients (default: 0, ordinary least squares)")
    args = parser.parse_args()
    source, destination = args.input_root.resolve(), args.output_root.resolve()
    if source == destination or source in destination.parents:
        parser.error("Output must be outside the input dataset")
    if args.components < 1 or not np.isfinite(args.ridge) or args.ridge < 0:
        parser.error("Components must be positive and ridge must be finite and nonnegative")
    selected = set(args.subject or [])
    scans = [p for p in sorted(source.glob("sub-*/ses-*/func"))
             if not selected or p.parent.parent.name in selected]
    missing = selected - {p.parent.parent.name for p in scans}
    if not scans or missing:
        parser.error(f"No matching inputs in {source}; unmatched subject IDs: {sorted(missing)}")

    completed = skipped = failed = 0
    for func in scans:
        relative = func.relative_to(source)
        output = destination / relative
        if output.exists():
            # Existing scientific results are never overwritten or silently reused
            # when the requested model differs.
            try:
                qc = json.loads((output / "regression_qc.json").read_text())
                if ((qc["aCompCor_components"], qc["ridge_penalty"]) != (args.components, args.ridge)
                        or len(data_files(output, "*_nuisance-regressed_bold.nii.gz")) != 1
                        or not (output / "nuisance_regressors.tsv").is_file()):
                    raise ValueError("existing output is incomplete or uses different settings")
            except (OSError, ValueError, KeyError) as error:
                failed += 1
                print(f"FAILED {relative}: {error}", flush=True)
                continue
            skipped += 1
            print(f"SKIPPED {relative}: complete output already exists", flush=True)
            continue
        try:
            run_scan(func, output, args.components, args.ridge)
        except (OSError, ValueError, np.linalg.LinAlgError, nib.filebasedimages.ImageFileError) as error:
            failed += 1
            print(f"FAILED {relative}: {error}", flush=True)
        else:
            completed += 1
            print(f"OK {relative}", flush=True)
    print(f"Completed {completed}; skipped {skipped}; failed {failed}", flush=True)
    return int(failed > 0)


if __name__ == "__main__":
    raise SystemExit(main())
