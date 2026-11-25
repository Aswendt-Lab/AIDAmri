"""Read fMRI metadata and select a repetition time in seconds."""

import argparse
import json
import math
from pathlib import Path


DEFAULT_TR_SECONDS = 1.42
LEGACY_TIMING_KEYS = {"ObjOrderList", "n_slices", "costum_timings"}


def positive_tr(value):
    """Validate a finite, positive TR; also suitable as an argparse type."""
    try:
        if isinstance(value, bool):
            raise ValueError
        tr = float(value)
        if not math.isfinite(tr) or tr <= 0:
            raise ValueError
    except (TypeError, ValueError, OverflowError):
        raise argparse.ArgumentTypeError("TR must be a finite number greater than zero (seconds)") from None
    return tr


def read_fmri_metadata(input_file):
    """Read the adjacent JSON for .nii/.nii.gz, returning {} on failure."""
    path = Path(input_file)
    if path.name.endswith(".nii.gz"):
        path = path.with_suffix("")
    sidecar = path.with_suffix(".json")
    try:
        with sidecar.open(encoding="utf-8") as infile:
            metadata = json.load(infile)
        if not isinstance(metadata, dict):
            raise ValueError("expected a JSON object")
    except (OSError, ValueError) as exc:
        print(f"Warning: Cannot use metadata JSON '{sidecar}': {exc}", flush=True)
        return {}
    return metadata


def resolve_tr(metadata, cli_tr=None):
    """Select and log TR: explicit CLI value, JSON, then 1.42 seconds.

    Legacy AIDAmri converters write milliseconds alongside the three custom
    timing keys. Other sidecars use BIDS RepetitionTime in seconds.
    """
    if cli_tr is not None:
        tr, source = positive_tr(cli_tr), "cli-option"
    else:
        tr, source = DEFAULT_TR_SECONDS, "fallback"
        if "RepetitionTime" in metadata:
            try:
                value = metadata["RepetitionTime"]
                if not isinstance(value, (int, float)) or isinstance(value, bool):
                    raise argparse.ArgumentTypeError("RepetitionTime must be a JSON number")
                json_tr = positive_tr(value)
                if LEGACY_TIMING_KEYS.issubset(metadata):
                    json_tr = positive_tr(json_tr / 1000.0)
                    print("Legacy AIDAmri JSON: converting RepetitionTime from milliseconds to seconds.", flush=True)
                tr, source = json_tr, "JSON"
            except argparse.ArgumentTypeError as exc:
                print(f"Warning: Invalid JSON RepetitionTime: {exc}; using fallback.", flush=True)
        else:
            print("Warning: No JSON RepetitionTime available; using fallback.", flush=True)

    print(f"Using TR = {tr} s (source: {source})", flush=True)
    return tr
