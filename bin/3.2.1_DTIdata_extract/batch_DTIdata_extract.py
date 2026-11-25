"""Run DTIdata_extract.py across proc_folder/subject/session/dwi folders.

Match the ROI name passed with -r immediately before .nii or .nii.gz,
ignoring case (e.g. annosplit does not match AnnoSplit_parental).
For each session, run the extractor sequentially on the selected maps
in DSI_studio. Select FA, MD, RD or AD with -m/--metric. The default
selection, all, processes all four metrics and legacy .fa0 maps.
Forward an optional -t translation table without choosing one automatically.

With --exclude-stroke-mask, search directly in each dwi folder for files ending
in Stroke_mask.nii or Stroke_mask.nii.gz, ignoring case and excluding subfolders.
Check all sessions before starting extractions; multiple matching masks or a
directory read error abort the batch. Forward a unique mask path to the extractor.
Without a mask, warn once per session and extract normally without the mask option.
Individual CSVs gain _exclude-stroke only when a mask is used. Cohort filenames
are unchanged and may combine sessions with and without stroke exclusion.
Retain missing and fully excluded ROI medians as empty cells in cohort tables.

The extractor creates individual CSV files beside the parameter maps.
Optionally pass -o/--output /OUTPUT/PATH/NAME (requires -t) to collect successful
exports from this run into /OUTPUT/PATH/NAME_<metric>.csv, one per metric.
An optional trailing .csv is removed before adding the metric suffix.
Relative output paths use the invocation's working directory; ~ is expanded.
Missing output directories are created automatically.
With -g/--group CSV (requires -o), insert a group column after session.
Group CSV headers name the groups; cells below contain subject folder names.
Ignore empty cells, match subjects exactly and leave missing groups empty.
Rows contain subject, session and median values in acronym columns ordered
as first encountered in the individual CSVs, following translation-table order;
append newly encountered acronyms at the end. Missing regions are empty cells,
including columns and sessions with all medians empty. Multiple images of one metric
in a session are ambiguous: report an error and omit that cohort row while
still running individual extractions. FA and legacy FA0 count as one metric.
Replace cohort tables atomically; leave existing tables unchanged with a
warning when no valid rows are available for a metric.
Cohort CSVs use normal new-file permissions, respecting umask and inherited permissions.
Show warnings and errors in the terminal; the extractor appends to the session's
analysis log, creating it if needed. Leave existing process logs unchanged.
Track dwi/extraction outputs around each extractor invocation and finalize
even after child-process failure. Pass the tracked folder only to that child
so the extractor does not start a second tracker for the same folder.
Suppress normal extractor output such as Done messages during batch runs.
Write one extractor header per session and batch invocation before the first
metric, inside extraction tracking; suppress the headers in all children.
Suppress the child's Git block and append Git information once per attempted
session to its existing analysis log after all metrics, including failures.
Fully skipped sessions do not receive a Git block.
List failed and fully or partially skipped sessions in the final summary.
The skipped counter counts skipped entries, which may be whole sessions
or individual missing metrics, rather than unique skipped sessions.
Continue after individual failures. Exit with status 1 if any extraction
or session fails, or no extraction succeeds; otherwise exit with status 0.

Usage:
    python batch_DTIdata_extract.py -i PROC_FOLDER -r annosplit_parental [-t TXT] [-m fa] [-o /OUTPUT/PATH/NAME] [-g GROUPS.csv] [--exclude-stroke-mask]

Created on 08.04.2019 by Niklas Pallast and Markus Aswendt.
latest Update: Julian Cramer
"""

import argparse
import csv
from dataclasses import dataclass
import math
import os
from pathlib import Path
import subprocess
import sys
from typing import NamedTuple
import uuid

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.artifact_manifest import manifest_filename, start_output_tracking
from common.dti_csv import getOutfile
from common.script_logging import append_git_information, build_script_log_path, format_script_log_header


METRICS = ("fa", "md", "rd", "ad")
PARENT_TRACKING_ENV = "AIDAMRI_DTI_EXTRACTION_TRACKED_BY_PARENT"
PARENT_GIT_LOG_ENV = "AIDAMRI_DTI_GIT_LOGGED_BY_PARENT"
PARENT_HEADER_LOG_ENV = "AIDAMRI_DTI_HEADER_LOGGED_BY_PARENT"


@dataclass
class _SessionLog:
    path: Path
    header_written: bool = False

    def write_header(self, context):
        if self.header_written:
            return True
        try:
            with self.path.open("a", encoding="utf-8") as log:
                log.write(f"\n{format_script_log_header('DTIdata_extract.py')}\n")
        except (OSError, UnicodeError) as error:
            print(f"Error: [{context}] Cannot write session header to '{self.path}': {error}", file=sys.stderr)
            return False
        self.header_written = True
        return True


class SessionResult(NamedTuple):
    succeeded: int
    skipped: int
    failed: int
    skip_reason: str = ""


def _nifti_stem(path):
    for suffix in (".nii.gz", ".nii"):
        if path.name.lower().endswith(suffix):
            return path.name[:-len(suffix)]
    return None


def findData(path):
    """Return sorted session DWI directories without assuming naming prefixes."""
    return sorted(folder for folder in Path(path).glob("*/*/dwi") if folder.is_dir())


def _find_rois(dwi, roi_name):
    return sorted(
        path for path in dwi.iterdir()
        if path.is_file() and (stem := _nifti_stem(path)) is not None
        and stem.lower().endswith(roi_name.lower())
    )


def _find_metrics(folder, selected_metric="all"):
    selected = METRICS if selected_metric == "all" else (selected_metric,)
    metrics = {metric: [] for metric in selected}
    for path in sorted(folder.iterdir()):
        stem = _nifti_stem(path)
        if not path.is_file() or stem is None:
            continue
        metric = stem.rsplit(".", 1)[-1].lower()
        if metric == "fa0" and selected_metric == "all":
            metric = "fa"
        if metric in metrics:
            metrics[metric].append(path)
    return metrics


def _cohort_output_base(name):
    """Return an absolute output base path without an optional CSV extension."""
    name = name.strip()
    if not name or "\0" in name or name.endswith(os.sep) or name.rsplit(os.sep, 1)[-1] in {".", ".."}:
        raise ValueError("--output must contain an output base filename, optionally preceded by a directory.")
    path = Path(name).expanduser()
    if path.is_dir():
        raise ValueError(f"--output must include a base filename, not an existing directory: '{path}'.")
    stem = path.name[:-4] if path.name.lower().endswith(".csv") else path.name
    if not stem.strip():
        raise ValueError("--output must contain a filename before the .csv extension.")
    return path.with_name(stem).resolve()


def _read_groups(path):
    """Map subject folder names to the group named by their CSV column header."""
    with path.open(encoding="utf-8-sig", newline="") as stream:
        preview = stream.read(8192)
        stream.seek(0)
        try:
            dialect = csv.Sniffer().sniff(preview, delimiters=",;\t")
        except csv.Error:
            # Unequal group sizes can make sniffing the full CSV inconclusive.
            try:
                dialect = csv.Sniffer().sniff(preview.splitlines()[0] if preview else "", delimiters=",;\t")
            except csv.Error:
                dialect = csv.excel
        reader = csv.reader(stream, dialect, strict=True)
        headers = [value.strip() for value in next(reader, [])]
        if not headers or any(not value for value in headers):
            raise ValueError("Group CSV must start with nonempty group column names.")
        if len(headers) != len(set(headers)):
            raise ValueError("Group CSV has duplicate group column names.")
        groups = {}
        for row in reader:
            if len(row) > len(headers):
                raise ValueError(f"Group CSV line {reader.line_num} has too many columns.")
            for group, subject in zip(headers, row):
                subject = subject.strip()
                if not subject:
                    continue
                if subject in groups:
                    raise ValueError(f"Subject '{subject}' occurs more than once in the group CSV.")
                groups[subject] = group
    return groups


def _read_region_values(path, include_group=False):
    """Retain median text, including blanks for missing or fully excluded ROIs."""
    values = {}
    roi_ids = set()
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, strict=True)
        if next(reader, None) != ["roi_id", "acronym", "median"]:
            raise ValueError("Expected CSV columns roi_id,acronym,median.")
        for row in reader:
            if len(row) != 3:
                raise ValueError(f"Expected three columns at CSV line {reader.line_num}.")
            roi_id = int(row[0])
            acronym, median = row[1:]
            if not acronym.strip() or acronym in {"subject", "session"} or (include_group and acronym == "group"):
                raise ValueError(f"Empty or reserved acronym '{acronym}' at CSV line {reader.line_num}.")
            if acronym in values or roi_id in roi_ids:
                raise ValueError(f"Ambiguous acronym or ROI ID at CSV line {reader.line_num}: '{acronym}', {roi_id}.")
            if median != "" and not math.isfinite(float(median)):
                raise ValueError(f"Nonfinite median for '{acronym}' at CSV line {reader.line_num}.")
            roi_ids.add(roi_id)
            values[acronym] = median
    return values


def _write_cohort_csv(path, rows, groups=None):
    """Replace one cohort table only after the entire new CSV has been written."""
    acronyms = list(dict.fromkeys(
        acronym
        for values in rows.values()
        for acronym in values
    ))
    temporary = None
    try:
        # Exclusive creation respects umask and inherited permissions.
        while True:
            candidate = path.parent / f".{path.name}.{uuid.uuid4().hex}.tmp"
            try:
                stream = candidate.open("x", encoding="utf-8", newline="")
            except FileExistsError:
                continue
            temporary = candidate
            break
        with stream:
            writer = csv.writer(stream)
            metadata = ["subject", "session"]
            if groups is not None:
                metadata.append("group")
            writer.writerow([*metadata, *acronyms])
            for (subject, session), values in sorted(rows.items()):
                identifiers = [subject, session]
                if groups is not None:
                    identifiers.append(groups.get(subject, ""))
                writer.writerow([*identifiers, *(values.get(acronym, "") for acronym in acronyms)])
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError as error:
                print(f"Warning: Cannot remove temporary cohort CSV '{temporary}': {error}", file=sys.stderr)


def _write_cohort_tables(output_base, cohort_rows, groups=None):
    written = []
    failed = 0
    for metric, rows in cohort_rows.items():
        path = output_base.with_name(f"{output_base.name}_{metric}.csv")
        if not rows:
            print(f"Warning: No valid {metric.upper()} rows for '{path}'. "
                  "No cohort CSV written; any existing file is unchanged.", file=sys.stderr)
            continue
        try:
            _write_cohort_csv(path, rows, groups)
        except (OSError, ValueError, csv.Error) as error:
            print(f"Error: Cannot write cohort CSV '{path}': {error}", file=sys.stderr)
            failed += 1
        else:
            written.append(str(path))
    return written, failed


def _print_summary(succeeded, skipped, failed, failed_sessions=None, reason="", skipped_sessions=None):
    message = f"Summary: {succeeded} succeeded, {skipped} skipped, {failed} failed."
    print(f"{message} {reason}" if reason else message, flush=True)
    if failed_sessions:
        print("Failed sessions:")
        for session, note in sorted(failed_sessions.items()):
            print(f"  {session}" + (f" ({note})" if note else ""))
    if skipped_sessions:
        print("Skipped sessions:")
        for session, note in sorted(skipped_sessions.items()):
            print(f"  {session} ({note})")
    return 1 if failed or not succeeded else 0


def _run_tracked_extraction(command, dwi, context, image, session_log):
    child_env = os.environ.copy()
    # Git logging belongs to the batch even without a recognizable BIDS manifest.
    child_env[PARENT_GIT_LOG_ENV] = str(dwi.resolve())
    child_env[PARENT_HEADER_LOG_ENV] = str(dwi.resolve())
    # An inherited marker must not suppress tracking without a matching owner.
    child_env.pop(PARENT_TRACKING_ENV, None)
    tracker = None
    try:
        manifest_filename(dwi, "dwi")
    except ValueError:
        # Keep support for arbitrary folder names; the extractor warns as before.
        pass
    except OSError as error:
        print(f"Error: [{context}] Cannot start manifest tracking for '{image.name}': {error}", file=sys.stderr)
        return False
    else:
        try:
            tracker = start_output_tracking(dwi, "dwi", "extraction")
        except (OSError, ValueError) as error:
            print(f"Error: [{context}] Cannot start manifest tracking for '{image.name}': {error}", file=sys.stderr)
            return False
        child_env[PARENT_TRACKING_ENV] = str(dwi.resolve())

    succeeded = False
    try:
        # Snapshot first so a log created by the batch is tracked as a new output.
        if not session_log.write_header(context):
            return False
        # The extractor owns its log. Keep stderr visible and suppress normal stdout.
        result = subprocess.run(command, cwd=dwi, check=False, stdout=subprocess.DEVNULL,
                                stderr=None, env=child_env)
        if result.returncode != 0:
            print(f"Error: [{context}] Extraction failed for '{image.name}' "
                  f"(exit code {result.returncode}).", file=sys.stderr)
        else:
            succeeded = True
    except OSError as error:
        print(f"Error: [{context}] Could not run extractor for '{image.name}': {error}", file=sys.stderr)
    finally:
        if tracker is not None:
            try:
                tracker.finalize()
            except (OSError, ValueError) as error:
                print(f"Error: [{context}] Cannot finalize extraction manifest for '{image.name}': {error}",
                      file=sys.stderr)
                succeeded = False
    return succeeded


def _append_session_git_information(dwi, context):
    log_path = Path(build_script_log_path(dwi, "analysis.log"))
    try:
        # Append only to an existing session log.
        if not log_path.exists():
            return True
        with log_path.open("a", encoding="utf-8") as log:
            append_git_information(log, output_dir=dwi)
    except (OSError, UnicodeError) as error:
        print(f"Error: [{context}] Cannot append Git information to '{log_path}': {error}", file=sys.stderr)
        return False
    return True


def _run_session(dwi, context, roi_name, metric_selection, translator, extractor, cohort_rows=None,
                 include_group=False, stroke_mask=None):
    succeeded = skipped = failed = 0
    missing_metrics = []
    try:
        rois = _find_rois(dwi, roi_name)
        if not rois:
            print(f"Warning: [{context}] No ROI matching '*{roi_name}.nii[.gz]'. Skipping session.",
                  file=sys.stderr)
            return SessionResult(0, 1, 0, f"fully skipped: ROI matching '*{roi_name}.nii[.gz]' missing")
        if len(rois) > 1:
            print(f"Error: [{context}] Multiple matching ROIs: {', '.join(p.name for p in rois)}. "
                  "Skipping session.", file=sys.stderr)
            return SessionResult(0, 0, 1)
        metrics_folder = dwi / "DSI_studio"
        if not metrics_folder.is_dir():
            print(f"Warning: [{context}] DSI_studio folder missing. Skipping session.", file=sys.stderr)
            return SessionResult(0, 1, 0, "fully skipped: DSI_studio folder missing")
        metrics = _find_metrics(metrics_folder, metric_selection)
    except OSError as error:
        print(f"Error: [{context}] Cannot read session inputs: {error}", file=sys.stderr)
        return SessionResult(0, 0, 1)

    # Different NIfTI encodings of one image would target the same CSV.
    by_stem = {}
    for images in metrics.values():
        for image in images:
            by_stem.setdefault(_nifti_stem(image), []).append(image)
    collisions = set()
    for stem, images in by_stem.items():
        if len(images) > 1:
            print(f"Error: [{context}] Colliding CSV outputs for '{stem}': "
                  f"{', '.join(p.name for p in images)}. Skipping these maps.", file=sys.stderr)
            collisions.update(images)
            failed += len(images)

    ambiguous_metrics = set()
    if cohort_rows is not None:
        for metric, images in metrics.items():
            if len(images) > 1:
                ambiguous_metrics.add(metric)
                print(f"Error: [{context}] Multiple {metric.upper()} images: {', '.join(p.name for p in images)}. "
                      "Omitting this subject/session/metric from the cohort CSV.", file=sys.stderr)
                failed += 1

    session_log = _SessionLog(Path(build_script_log_path(dwi, "analysis.log")))
    try:
        for metric, images in metrics.items():
            if not images:
                print(f"Warning: [{context}] No {metric.upper()} map found. Skipping parameter.", file=sys.stderr)
                skipped += 1
                missing_metrics.append(metric.upper())
            for image in images:
                if image in collisions:
                    continue
                command = [sys.executable, str(extractor), "-i", str(image), "-r", str(rois[0])]
                if translator is not None:
                    command.extend(("-t", str(translator)))
                if stroke_mask is not None:
                    command.extend(("--exclude-stroke-mask", str(stroke_mask)))
                if _run_tracked_extraction(command, dwi, context, image, session_log):
                    succeeded += 1
                    if cohort_rows is not None and metric not in ambiguous_metrics:
                        csv_path = Path(getOutfile(rois[0], image, exclude_stroke=stroke_mask is not None))
                        try:
                            values = _read_region_values(csv_path, include_group=include_group)
                        except (OSError, ValueError, csv.Error) as error:
                            print(f"Error: [{context}] Cannot collect '{csv_path}' for the cohort CSV: {error}",
                                  file=sys.stderr)
                            failed += 1
                        else:
                            if values:
                                key = (dwi.parent.parent.name, dwi.parent.name)
                                cohort_rows[metric][key] = values
                            else:
                                print(f"Warning: [{context}] No region values in '{csv_path}'. "
                                      "Omitting this cohort row.", file=sys.stderr)
                else:
                    failed += 1
    finally:
        if session_log.header_written and not _append_session_git_information(dwi, context):
            failed += 1
    skip_reason = ""
    if missing_metrics:
        status = "fully" if succeeded == 0 and failed == 0 else "partially"
        skip_reason = f"{status} skipped: {', '.join(missing_metrics)} missing"
    return SessionResult(succeeded, skipped, failed, skip_reason)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Batch script controlling DTIdata_extract.py across subject/session/dwi folders."
                                                 "Summarized infos in cohort csvs.")
    parser.add_argument("-i", "--input", required=True,
                        help="Cohort root containing subject/session/dwi folders")
    parser.add_argument("-r", "--roi_name", required=True,
                        help="ROI filename suffix before .nii[.gz], e.g. annosplit or annosplit_parental")
    parser.add_argument("-t", "--translatorTXT",
                        help="Translation TXT defining ROI rows and their order; regions absent from an "
                             "atlas have empty medians. Passed to DTIdata_extract.py; required with --output")
    parser.add_argument("-m", "--metric", choices=(*METRICS, "all"), default="all", type=str.lower,
                        help="Metric suffix before .nii[.gz] (default: all; includes legacy .fa0 maps)")
    parser.add_argument("-o", "--output", metavar="PATH/NAME",
                        help="Output base path: write PATH/NAME_<metric>.csv; "
                             "requires -t")
    parser.add_argument("-g", "--group", metavar="CSV",
                        help="Group CSV: headers are group names, cells are subject folder names; requires -o")
    parser.add_argument(
        "--exclude-stroke-mask", action="store_true",
        help="Search only directly in each dwi folder for *Stroke_mask.nii or .nii.gz "
             "(case insensitive). Use a unique mask for ROI medians; warn and extract normally "
             "if absent. Multiple matches abort the batch before any extraction. "
             "Cohort CSV filenames are unchanged.",
    )
    args = parser.parse_args(argv)

    root = Path(args.input).expanduser().resolve()
    translator = Path(args.translatorTXT).expanduser().resolve() if args.translatorTXT is not None else None
    group_file = Path(args.group).expanduser().resolve() if args.group is not None else None
    extractor = Path(__file__).resolve().with_name("DTIdata_extract.py")
    roi_name = args.roi_name.strip()
    if not root.is_dir():
        return _print_summary(0, 0, 1, reason=f"Error: Proc folder does not exist or is not a directory: '{root}'.")
    if not roi_name:
        return _print_summary(0, 0, 1, reason="Error: ROI name must not be empty.")
    if group_file is not None and args.output is None:
        return _print_summary(0, 0, 1, reason="Error: --group requires -o/--output for cohort CSVs.")
    output_base = None
    cohort_rows = None
    if args.output is not None:
        if translator is None:
            return _print_summary(0, 0, 1, reason="Error: --output requires -t/--translatorTXT for acronym columns.")
        try:
            output_base = _cohort_output_base(args.output)
        except (OSError, ValueError, RuntimeError) as error:
            return _print_summary(0, 0, 1, reason=f"Error: {error}")
        selected = METRICS if args.metric == "all" else (args.metric,)
        cohort_rows = {metric: {} for metric in selected}
    for path in (extractor, translator, group_file):
        if path is not None and not path.is_file():
            return _print_summary(0, 0, 1, reason=f"Error: Required file does not exist: '{path}'.")
    groups = None
    if group_file is not None:
        try:
            groups = _read_groups(group_file)
        except (OSError, ValueError, csv.Error) as error:
            return _print_summary(0, 0, 1, reason=f"Error: Cannot read group CSV '{group_file}': {error}")

    succeeded = skipped = failed = 0
    failed_sessions = {}
    skipped_sessions = {}
    try:
        sessions = findData(root)
    except OSError as error:
        return _print_summary(0, 0, 1, reason=f"Error: Cannot search '{root}': {error}")
    if not sessions:
        return _print_summary(0, 0, 0, reason=f"No subject/session/dwi folders found in '{root}'.")

    stroke_masks = {}
    if args.exclude_stroke_mask:
        # Check every session before starting children or creating output folders.
        for dwi in sessions:
            context = str(dwi.relative_to(root))
            try:
                masks = _find_rois(dwi, "Stroke_mask")
            except OSError as error:
                print(f"Error: [{context}] Cannot search for stroke masks: {error}. Aborting batch.",
                      file=sys.stderr)
                return _print_summary(0, 0, 1)
            if len(masks) > 1:
                print(f"Error: [{context}] Multiple matching stroke masks: "
                      f"{', '.join(path.name for path in masks)}. Aborting batch.", file=sys.stderr)
                session = f"{dwi.parent.parent.name} / {dwi.parent.name}"
                return _print_summary(0, 0, 1, failed_sessions={session: "multiple stroke masks"})
            stroke_masks[dwi] = masks[0] if masks else None
            if not masks:
                print(f"Warning: [{context}] No stroke mask matching '*Stroke_mask.nii[.gz]' "
                      "directly in dwi. Calculating medians without stroke exclusion.", file=sys.stderr)

    if output_base is not None:
        try:
            for metric in cohort_rows:
                path = output_base.with_name(f"{output_base.name}_{metric}.csv")
                if path.is_dir():
                    raise ValueError(f"Cohort CSV target is a directory: '{path}'.")
            output_base.parent.mkdir(parents=True, exist_ok=True)
        except (OSError, ValueError) as error:
            return _print_summary(0, 0, 1, reason=f"Error: Cannot prepare cohort output '{output_base}': {error}")

    for dwi in sessions:
        context = str(dwi.relative_to(root))
        session = f"{dwi.parent.parent.name} / {dwi.parent.name}"
        result = _run_session(dwi, context, roi_name, args.metric, translator, extractor, cohort_rows,
                              include_group=groups is not None, stroke_mask=stroke_masks.get(dwi))
        succeeded += result.succeeded
        skipped += result.skipped
        failed += result.failed
        if result.failed:
            failed_sessions.setdefault(session, "")
        if result.skip_reason:
            skipped_sessions[session] = result.skip_reason

    reason = ""
    if cohort_rows is not None:
        written, write_failures = _write_cohort_tables(output_base, cohort_rows, groups)
        failed += write_failures
        if written:
            reason = "Cohort CSVs: " + ", ".join(written) + "."
    return _print_summary(succeeded, skipped, failed, failed_sessions, reason=reason, skipped_sessions=skipped_sessions)


if __name__ == "__main__":
    sys.exit(main())
