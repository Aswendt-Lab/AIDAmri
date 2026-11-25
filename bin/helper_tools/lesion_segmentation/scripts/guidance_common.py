"""Lightweight shared geometry and metrics for the portable guidance workflow."""
import csv
import numpy as np


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sphere_mask(shape, spacing, regions):
    result = np.zeros(shape, bool)
    for region in regions:
        center = np.array(region['point']); radius = region['radius_mm']
        widths = np.ceil(radius / np.asarray(spacing)).astype(int)
        lo = np.maximum(center - widths, 0); hi = np.minimum(center + widths + 1, shape)
        slices = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
        grid = np.ogrid[tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))]
        dist = sum(((g-c)*s)**2 for g, c, s in zip(grid, center, spacing))
        result[slices] |= dist <= radius**2 + 1e-10
    return result


def metric(pred, ref, voxel):
    predicted = int(pred.sum())
    if ref is None:
        return dict(dice=None,precision=None,recall=None,predicted_volume_mm3=float(predicted*voxel),reference_volume_mm3=None,absolute_volume_error_mm3=None,true_positive_voxels=None,false_positive_voxels=None,false_negative_voxels=None,containment_ratio=None,prediction_containment_ratio=None)
    tp = int(np.count_nonzero(pred & ref)); fp = int(np.count_nonzero(pred & ~ref)); fn = int(np.count_nonzero(~pred & ref)); reference = int(ref.sum())
    return dict(dice=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 1.,precision=tp/(tp+fp) if tp+fp else 0.,recall=tp/(tp+fn) if tp+fn else 0.,predicted_volume_mm3=float(predicted*voxel),reference_volume_mm3=float(reference*voxel),absolute_volume_error_mm3=float(abs(predicted-reference)*voxel),true_positive_voxels=tp,false_positive_voxels=fp,false_negative_voxels=fn,containment_ratio=tp/reference if reference else None,prediction_containment_ratio=tp/predicted if predicted else None)
