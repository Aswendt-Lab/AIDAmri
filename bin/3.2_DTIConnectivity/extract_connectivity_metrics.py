"""Export one DSI connectivity metric across a BIDS-like cohort to CSV.

Read existing MAT files only. Each output row represents one subject/session;
columns contain individual regions (t2r) or unique region pairs (r2r).
Sessions missing the selected MAT file are skipped with a warning.
Exported CSV files use the system's normal permissions for new files.
"""

import argparse
import csv
from dataclasses import dataclass
import difflib
from itertools import chain, combinations
from pathlib import Path
import re
import sys
import uuid
import os

import numpy as np
from scipy.io import loadmat, whosmat


METADATA = ("subject", "session", "group")


@dataclass(frozen=True)
class Sample:
    subject: str
    session: str
    path: Path
    regions: tuple
    variable: str


@dataclass(frozen=True)
class ExportSummary:
    sessions: int
    subjects: int
    columns: int
    output: Path


def matlab_identifier(value):
    """Match the identifier conversion used by dsi_tools.py."""
    name = re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_") or "variable"
    if not name[0].isalpha():
        name = "v_" + name
    return name[:63]


def read_groups(path):
    """Column headers are group names; cells below them contain subject IDs."""
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        preview = stream.read(8192)
        stream.seek(0)
        try:
            dialect = csv.Sniffer().sniff(preview, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.reader(stream, dialect)
        headers = next(reader, [])
        headers = [value.strip() for value in headers]
        if not headers or any(not value for value in headers):
            raise ValueError("Group CSV must start with nonempty group column names.")
        if len(headers) != len(set(headers)):
            raise ValueError("Group CSV has duplicate group column names.")
        groups = {}
        for line, row in enumerate(reader, 2):
            if len(row) > len(headers):
                raise ValueError(f"Group CSV line {line} has too many columns.")
            for group, subject in zip(headers, row):
                subject = subject.strip()
                if not subject:
                    continue
                if subject in groups:
                    raise ValueError(f"Subject {subject!r} occurs more than once in the group CSV.")
                groups[subject] = group
        if not groups:
            raise ValueError("Group CSV contains no subjects.")
        return groups


def _text(value):
    array = np.asarray(value)
    if array.dtype.kind in "ui":
        return array.astype(np.uint8).ravel().tobytes().decode("utf-8")
    if array.dtype.kind == "U":
        return "".join(array.ravel())
    if array.dtype.kind == "S":
        return b"".join(array.ravel()).decode("utf-8")
    raise ValueError(f"Unsupported region-name type: {array.dtype}")


def read_regions(data):
    if "region_names" in data:
        regions = tuple(_text(value).rstrip("\x00") for value in data["region_names"].ravel())
    elif "name" in data:
        regions = tuple(_text(data["name"]).rstrip("\x00\r\n").splitlines())
    else:
        raise ValueError("MAT file contains neither 'name' nor 'region_names'.")
    if not regions or any(not region.strip() for region in regions):
        raise ValueError("MAT file contains empty region names.")
    if len(regions) != len(set(regions)):
        raise ValueError("MAT file contains duplicate region names.")
    return regions


def select_metric(variables, metric, level):
    requested = matlab_identifier(metric)
    available = [(name, shape) for name, shape, _ in variables
                 if matlab_identifier(name).endswith("_" + level)]
    matches = [(name, shape) for name, shape in available
               if matlab_identifier(name) == requested + "_" + level]
    if len(matches) > 1:
        raise ValueError(f"Metric {metric!r} is ambiguous: {[name for name, _ in matches]}")
    if not matches:
        choices = sorted({matlab_identifier(name)[:-(len(level) + 1)] for name, _ in available})
        raise ValueError(f"Metric {metric!r} is unavailable for {level}. Available metrics: "
                         + ", ".join(choices))
    return matches[0]


def discover_samples(root, atlas, connectivity, level, metric):
    root = Path(root)
    if not root.is_dir():
        raise ValueError(f"Input folder does not exist: {root}")
    directories = sorted(root.glob("sub-*/ses-*/dwi"))
    if not directories:
        raise ValueError(f"No sub-*/ses-*/dwi folders found in {root}")
    marker = "AnnoSplit_parental" if atlas == "parental" else "AnnoSplit"
    pattern = re.compile(rf"{marker}(?:_resampled)?\.all\.{connectivity}\.connectivity\.mat$")
    samples = []
    for directory in directories:
        candidates = sorted(path for path in (directory / "connectivity").glob("*.mat")
                            if pattern.search(path.name))
        if not candidates:
            print(f"WARNING: Skipping {directory.parent.parent.name}/{directory.parent.name}: "
                  f"no {atlas}/{connectivity} MAT file found in {directory / 'connectivity'}.",
                  file=sys.stderr)
            continue
        if len(candidates) != 1:
            raise ValueError(f"Expected one {atlas}/{connectivity} MAT file in {directory}; "
                             f"found {len(candidates)}.")
        path = candidates[0]
        try:
            regions = read_regions(loadmat(path, variable_names=["name", "region_names"]))
            variable, shape = select_metric(whosmat(path), metric, level)
            n = len(regions)
            expected = [(n, 1), (1, n)] if level == "t2r" else [(n, n)]
            if shape not in expected:
                raise ValueError(f"{variable!r} has shape {shape}; expected {expected}.")
        except (ValueError, OSError) as error:
            raise ValueError(f"{path}: {error}") from error
        samples.append(Sample(directory.parent.parent.name, directory.parent.name,
                              path, regions, variable))
    if not samples:
        raise ValueError(f"No sessions with a matching {atlas}/{connectivity} MAT file found in {root}.")
    return samples


def resolve_region(query, regions):
    """Accept only exact region names stored in the selected MAT files."""
    if query in regions:
        return query
    suggestions = difflib.get_close_matches(query, regions, n=5)
    hint = f" Similar names: {', '.join(suggestions)}." if suggestions else ""
    raise ValueError(f"Region {query!r} is absent from the selected atlas results. "
                     f"Enter a full region name exactly as stored in the MAT file.{hint}")


def read_values(sample, level, metric):
    values = loadmat(sample.path, variable_names=[sample.variable])[sample.variable]
    if values.dtype.kind not in "iuf":
        raise ValueError(f"{sample.path}: metric must contain real numeric values.")
    if level == "r2r" and not np.array_equal(values, values.T, equal_nan=True):
        raise ValueError(f"{sample.path}: r2r matrix is asymmetric; unique unordered pairs "
                         "would discard different values.")
    if matlab_identifier(metric) == "number_of_tracts":
        if not np.all(np.isfinite(values) & (values >= 0) & (values == np.floor(values))):
            raise ValueError(f"{sample.path}: track counts must be finite nonnegative integers.")
    return values.ravel() if level == "t2r" else values


def _number(value):
    return format(float(value), ".17g")


def row_values(sample, values, regions, level, pair):
    indices = {name: index for index, name in enumerate(sample.regions)}
    if pair:
        left, right = pair
        yield _number(values[indices[left], indices[right]]) if left in indices and right in indices else ""
    elif level == "t2r":
        for name in regions:
            yield _number(values[indices[name]]) if name in indices else ""
    else:
        for i, left in enumerate(regions):
            left_index = indices.get(left)
            for j in range(i + 1, len(regions)):
                right_index = indices.get(regions[j])
                yield (_number(values[left_index, right_index])
                       if left_index is not None and right_index is not None else "")


def export_cohort(root, group_csv, output, atlas, connectivity, level,
                  metric="number_of_tracts", region1=None, region2=None):
    if atlas not in ("parental", "detailed") or connectivity not in ("pass", "end") or level not in ("t2r", "r2r"):
        raise ValueError("Invalid atlas, connectivity or level selection.")
    if bool(region1) != bool(region2):
        raise ValueError("Provide both --region1 and --region2 together.")
    if region1 and level != "r2r":
        raise ValueError("--region1 and --region2 are available only for --level r2r.")
    if not metric.strip():
        raise ValueError("Metric must not be empty.")
    output = Path(output).resolve()
    if output.suffix.lower() != ".csv":
        raise ValueError("Output path must end in .csv.")
    if group_csv is not None and output == Path(group_csv).resolve():
        raise ValueError("Output must not overwrite the group CSV.")

    metadata = METADATA if group_csv is not None else METADATA[:2]
    groups = read_groups(group_csv) if group_csv is not None else {}
    samples = discover_samples(root, atlas, connectivity, level, metric)
    missing_groups = sorted({sample.subject for sample in samples} - groups.keys())
    if group_csv is not None and missing_groups:
        raise ValueError("Subjects missing from group CSV: " + ", ".join(missing_groups))
    regions = sorted({name for sample in samples for name in sample.regions})
    pair = None
    if region1:
        pair = (resolve_region(region1, regions), resolve_region(region2, regions))
        if pair[0] == pair[1]:
            raise ValueError("Select two different regions; diagonal values are excluded.")
        pair = tuple(sorted(pair))
    if level == "r2r" and any("_x_" in name for name in regions):
        raise ValueError("A region name contains the reserved pair separator '_x_'.")
    if level == "t2r" and set(regions).intersection(metadata):
        raise ValueError("A region name conflicts with a subject/session/group column.")

    if pair:
        count = 1
        columns = iter(["_x_".join(pair)])
    elif level == "t2r":
        count = len(regions)
        columns = iter(regions)
    else:
        count = len(regions) * (len(regions) - 1) // 2
        columns = ("_x_".join(pair) for pair in combinations(regions, 2))

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        # Exclusive creation respects umask and inherited permissions, while
        # keeping incomplete exports separate from the destination.
        while True:
            candidate = output.parent / f".connectivity_{uuid.uuid4().hex}.csv"
            try:
                stream = candidate.open("x", encoding="utf-8", newline="")
            except FileExistsError:
                continue
            temporary_path = candidate
            break
        with stream:
            writer = csv.writer(stream)
            writer.writerow(chain(metadata, columns))
            for sample in samples:
                values = read_values(sample, level, metric)
                identifiers = (sample.subject, sample.session)
                if group_csv is not None:
                    identifiers += (groups[sample.subject],)
                writer.writerow(chain(identifiers,
                                      row_values(sample, values, regions, level, pair)))
        os.replace(temporary_path, output)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()
    return ExportSummary(len(samples), len({sample.subject for sample in samples}),
                         len(metadata) + count, output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, required=True, help="BIDS-like proc folder")
    parser.add_argument("-a", "--atlas", choices=("parental", "detailed"), required=True)
    parser.add_argument("-c", "--connectivity", choices=("pass", "end"), required=True)
    parser.add_argument("-l", "--level", choices=("t2r", "r2r"), required=True,
                        help="One value per region or per unique region pair")
    parser.add_argument("-m", "--metric", default="number_of_tracts",
                        help="Metric without t2r/r2r suffix (default: number_of_tracts)")
    parser.add_argument("-g", "--group", type=Path,
                        help="Optional CSV: group names in the header, subject IDs below; "
                             "omit to exclude the group column")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Destination CSV")
    parser.add_argument("-r1", "--region1", help="Exact first r2r region name from the MAT file, e.g. L_Primary_motor_area")
    parser.add_argument("-r2", "--region2", help="Exact second r2r region name from the MAT file, e.g. R_Primary_motor_area")
    args = parser.parse_args(argv)
    if bool(args.region1) != bool(args.region2):
        parser.error("Provide both --region1 and --region2 together.")
    if args.region1 and args.level != "r2r":
        parser.error("--region1 and --region2 are available only for --level r2r.")
    try:
        summary = export_cohort(args.input, args.group, args.output, args.atlas,
                                args.connectivity, args.level, args.metric,
                                args.region1, args.region2)
    except (ValueError, OSError, NotImplementedError) as error:
        parser.error(str(error))
    print(f"Saved {summary.sessions} sessions from {summary.subjects} subjects, "
          f"{summary.columns} columns ({args.atlas}/{args.connectivity}/{args.level}, "
          f"metric={args.metric}) to {summary.output}")


if __name__ == "__main__":
    main()
