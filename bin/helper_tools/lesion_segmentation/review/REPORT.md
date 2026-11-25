# Saved human-guided review — 2026-11-25

User assessment: **working perfectly**. Audit found **12 saved runs / 11 distinct scans**, including one retry. Latest run per scan is retained for quantification; the full history remains preserved. All reviewed scans are P1; none are P7. 28 supplied scans have no saved guided run.

The 8 reviewed scans from the earlier strict-new subset are shown separately. The other three include two training cases and one prior validation image. This is post-validation manual refinement on a selected early subset, not an independent generalization estimate. Native direction remains primary; mirrored is reported separately without reference-based selection.

| Subset | Direction | Stage | n | Mean Dice | Containment | Prediction containment |
|---|---|---|---:|---:|---:|---:|
| all_reviewed | native | baseline | 11 | 0.454 | 93.6% | 30.5% |
| all_reviewed | native | guided_continuity | 11 | 0.836 | 93.2% | 77.4% |
| all_reviewed | mirrored | baseline | 11 | 0.454 | 93.6% | 30.5% |
| all_reviewed | mirrored | guided_continuity | 11 | 0.836 | 93.2% | 77.4% |
| strict_new_reviewed | native | baseline | 8 | 0.462 | 94.0% | 31.3% |
| strict_new_reviewed | native | guided_continuity | 8 | 0.823 | 93.6% | 75.4% |
| strict_new_reviewed | mirrored | baseline | 8 | 0.462 | 94.0% | 31.3% |
| strict_new_reviewed | mirrored | guided_continuity | 8 | 0.823 | 93.6% | 75.4% |

## Latest native results

| Scan | Dice baseline → guided | Guided containment | Guidance |
|---|---:|---:|---|
| GV_T3_3_1_1_2_20181121_082244 | 0.377 → 0.852 | 84.1% | 1 positive / 1 avoidance |
| GV_T3_3_2_1_2_20181121_091427 | 0.380 → 0.851 | 89.9% | 1 positive / 1 avoidance |
| GV_T3_3_3_1_2_20181121_102307 | 0.341 → 0.869 | 90.1% | 2 positive / 2 avoidance |
| GV_T3_4_1_1_3_20181121_125526 | 0.540 → 0.679 | 97.3% | 1 positive / 1 avoidance |
| GV_T3_4_2_1_2_20181121_133734 | 0.480 → 0.861 | 90.4% | 1 positive / 1 avoidance |
| GV_T3_4_3_1_2_20181121_143641 | 0.308 → 0.887 | 92.4% | 1 positive / 1 avoidance |
| GV_T3_8_2_1_2_20190206_092529 | 0.498 → 0.922 | 98.8% | 1 positive / 1 avoidance |
| GV_T3_8_4_1_2_20190206_104319 | 0.474 → 0.885 | 96.4% | 1 positive / 1 avoidance |
| GV_T3_9_1_1_2_20190206_112248 | 0.559 → 0.573 | 90.7% | 1 positive / 0 avoidance |
| GV_T3_9_3_1_2_20190206_130351 | 0.372 → 0.886 | 97.5% | 1 positive / 1 avoidance |
| SP_V1_2_1_1_1_20171013_144832 | 0.661 → 0.934 | 97.9% | 2 positive / 1 avoidance |

Containment = intersection/manual volume (recall); prediction containment = intersection/predicted volume (precision). See CSV for per-case volumes and summary medians/ranges.

One early retry produced an empty mask because its positive seed missed the prediction. The later saved run for that scan has a retained positive seed and one missed seed; its warning is preserved. No session was selected by a better Dice score, and no earlier failures were deleted. All latest native runs improve Dice over their baseline, but scores still vary: inspect the per-case report rather than interpreting the user assessment as perfect quantitative agreement.

All 24 saved masks across all 12 sessions were reproduced exactly from guidance, with geometry, 3D continuity, one-region-per-slice, avoidance exclusion and all 48 metric rows checked. Original images, reference arrays and frozen checkpoint remain unchanged. Guidance changes selection only; it cannot add tissue absent from the network output.

`selected_sessions.json`, `session_history.json`, `STATUS.json` and `source_checksums.json` freeze the audit. Originals remain in the neighbouring `manual_sessions/` folder.
