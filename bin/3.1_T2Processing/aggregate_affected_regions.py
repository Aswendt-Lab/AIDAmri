"""Combine existing affectedRegions CSVs across a BIDS-like proc folder.

Write one row per subject/session, with all regional stroke volumes first,
followed by all regional percentages in the same alphabetical region order.
Only regions appearing in at least one selected CSV become output columns.
Optionally translate Label_id values into hemisphere-prefixed acronyms.
"""

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import os
from pathlib import Path
import sys
import uuid


REGION = "Brain region"
VOLUME = "Region stroke volume [mm3]"
PERCENT = "Affected region percentage [%]"
TOTALS = ("Total stroke volume [mm3]", "Stroke percent of brain [%]")
REQUIRED = (REGION, VOLUME, PERCENT) + TOTALS


@dataclass(frozen=True)
class Sample:
    subject: str
    session: str
    path: Path
    regions: dict
    totals: tuple


@dataclass(frozen=True)
class ExportSummary:
    sessions: int
    subjects: int
    regions: int
    output: Path


def numeric_value(value, description):
    """Validate numbers without changing the source text or its precision."""
    try:
        number = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(f"{description}: expected a number, got {value!r}.") from error
    if not number.is_finite():
        raise ValueError(f"{description}: expected a finite number, got {value!r}.")
    return number


def read_translation(path):
    """Expand base ID/acronym pairs into left/right labels in memory only."""
    path = Path(path)
    base = {}
    acronyms = {}
    try:
        with path.open(encoding="utf-8-sig") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                fields = line.split(maxsplit=1)
                if len(fields) != 2:
                    raise ValueError(f"Line {line_number}: expected a label ID and acronym.")
                try:
                    label_id = int(fields[0])
                except ValueError as error:
                    raise ValueError(f"Line {line_number}: invalid integer label ID {fields[0]!r}.") from error
                acronym = fields[1].strip()
                if label_id in base:
                    raise ValueError(f"Line {line_number}: duplicate label ID {label_id}.")
                if acronym in acronyms:
                    raise ValueError(f"Line {line_number}: ambiguous acronym {acronym!r} for IDs "
                                     f"{acronyms[acronym]} and {label_id}.")
                base[label_id] = (acronym, line_number)
                acronyms[acronym] = label_id
        if not base:
            raise ValueError("Translation TXT contains no label ID/acronym pairs.")
        translation = {}
        for label_id, (acronym, line_number) in base.items():
            for expanded_id, prefix in ((label_id, "L_"), (label_id + 2000, "R_")):
                if expanded_id in translation:
                    raise ValueError(f"Line {line_number}: ID collision at {expanded_id} "
                                     "when expanding left/right labels.")
                translation[expanded_id] = prefix + acronym
    except (ValueError, UnicodeError) as error:
        raise ValueError(f"{path}: {error}") from error
    return translation


def read_sample(path, subject, session, translation=None):
    regions = {}
    totals = {}
    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream, strict=True)
            headers = reader.fieldnames or []
            required = REQUIRED + (("Label_id",) if translation is not None else ())
            missing = [name for name in required if name not in headers]
            if missing:
                raise ValueError("Missing required columns: " + ", ".join(missing))
            named_headers = [name for name in headers if name]
            if len(named_headers) != len(set(named_headers)):
                raise ValueError("Duplicate column names.")
            for row in reader:
                line = reader.line_num
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"Line {line}: column count differs from the header.")
                region = row[REGION]
                if not region.strip():
                    raise ValueError(f"Line {line}: empty brain region.")
                if translation is not None:
                    try:
                        label_id = int(row["Label_id"])
                    except ValueError as error:
                        raise ValueError(f"Line {line}: invalid integer Label_id {row['Label_id']!r}.") from error
                    if label_id not in translation:
                        raise ValueError(f"Line {line}: no translation for Label_id {label_id}.")
                    region = translation[label_id]
                    if region in regions:
                        raise ValueError(f"Line {line}: duplicate Label_id {label_id} ({region}).")
                if region in regions:
                    raise ValueError(f"Line {line}: duplicate brain region {region!r}.")
                values = tuple(row[column].strip() for column in (VOLUME, PERCENT))
                for column, value in zip((VOLUME, PERCENT), values):
                    numeric_value(value, f"Line {line}, {column}")
                regions[region] = values
                for column in TOTALS:
                    value = row[column].strip()
                    if not value:
                        continue
                    number = numeric_value(value, f"Line {line}, {column}")
                    if column in totals and number != Decimal(totals[column]):
                        raise ValueError(f"Line {line}: conflicting values for {column}.")
                    totals.setdefault(column, value)
        if not regions:
            raise ValueError("No data rows; total stroke values are unavailable.")
        missing = [name for name in TOTALS if name not in totals]
        if missing:
            raise ValueError("Missing total stroke values: " + ", ".join(missing))
    except (ValueError, csv.Error, UnicodeError) as error:
        raise ValueError(f"{path}: {error}") from error
    return Sample(subject, session, path, regions, tuple(totals[name] for name in TOTALS))


def discover_samples(root, suffix, translation=None):
    root = Path(root)
    if not root.is_dir():
        raise ValueError(f"Input folder does not exist: {root}")
    if not suffix or not suffix.lower().endswith(".csv") or "/" in suffix or "\\" in suffix:
        raise ValueError("Suffix must be a filename suffix ending in .csv, e.g. parental.csv.")
    samples = []
    for session in sorted(root.glob("sub-*/ses-*")):
        if not session.is_dir():
            continue
        directory = session / "anat" / "affected_Regions"
        matches = sorted(path for path in directory.glob("*")
                         if path.is_file() and "affectedregions" in path.name.lower()
                         and path.name.lower().endswith(suffix.lower()))
        if not matches:
            print(f"WARNING: Skipping {session.parent.name}/{session.name}: "
                  f"no affectedRegions CSV ending in {suffix!r} in {directory}.", file=sys.stderr)
            continue
        if len(matches) > 1:
            raise ValueError(f"Multiple matching CSVs in {directory}: "
                             + ", ".join(path.name for path in matches))
        samples.append(read_sample(matches[0], session.parent.name, session.name, translation))
    if not samples:
        raise ValueError(f"No matching affectedRegions CSVs ending in {suffix!r} found in {root}.")
    return samples


def read_groups(path):
    """Read subject IDs below group-name headers, including uneven columns."""
    path = Path(path)
    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            preview = stream.read(8192)
            stream.seek(0)
            try:
                dialect = csv.Sniffer().sniff(preview, delimiters=",;\t")
            except csv.Error:
                dialect = csv.excel
            reader = csv.reader(stream, dialect, strict=True)
            headers = [value.strip() for value in next(reader, [])]
            if not headers or any(not value for value in headers):
                raise ValueError("Group CSV must start with nonempty group column names.")
            if len(headers) != len(set(headers)):
                raise ValueError("Duplicate group column names.")
            groups = {}
            for row in reader:
                if len(row) > len(headers):
                    raise ValueError(f"Line {reader.line_num}: too many columns.")
                for group, subject in zip(headers, row):
                    subject = subject.strip()
                    if not subject:
                        continue
                    if subject in groups:
                        raise ValueError(f"Subject {subject!r} occurs more than once in the group CSV.")
                    groups[subject] = group
    except (ValueError, csv.Error, UnicodeError) as error:
        raise ValueError(f"{path}: {error}") from error
    return groups


def export_cohort(root, suffix, output, group_csv=None, translation_txt=None):
    translation = read_translation(translation_txt) if translation_txt is not None else None
    samples = discover_samples(root, suffix, translation)
    output = Path(output)
    sources = [sample.path for sample in samples]
    if group_csv is not None:
        sources.append(Path(group_csv))
    if translation_txt is not None:
        sources.append(Path(translation_txt))
    if any(output.resolve() == path.resolve() for path in sources):
        raise ValueError(f"Output must not overwrite an input file: {output}")
    groups = read_groups(group_csv) if group_csv is not None else {}
    subjects = {sample.subject for sample in samples}
    if group_csv is not None:
        for subject in sorted(subjects - groups.keys()):
            print(f"WARNING: No group assignment for {subject}; leaving group empty.", file=sys.stderr)
    regions = sorted({region for sample in samples for region in sample.regions})
    metadata = ["subject", "session"] + (["group"] if group_csv is not None else [])
    header = metadata + list(TOTALS)
    header += [f"{region} [mm3]" for region in regions]
    header += [f"{region} [%]" for region in regions]

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        # Exclusive creation respects normal file permissions and keeps failed
        # writes separate from any existing destination CSV.
        while True:
            candidate = output.parent / f".affected_regions_{uuid.uuid4().hex}.csv"
            try:
                stream = candidate.open("x", encoding="utf-8", newline="")
            except FileExistsError:
                continue
            temporary = candidate
            break
        with stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            for sample in samples:
                row = [sample.subject, sample.session]
                if group_csv is not None:
                    row.append(groups.get(sample.subject, ""))
                row.extend(sample.totals)
                row.extend(sample.regions.get(region, ("0", "0"))[0] for region in regions)
                row.extend(sample.regions.get(region, ("0", "0"))[1] for region in regions)
                writer.writerow(row)
        os.replace(temporary, output)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return ExportSummary(len(samples), len(subjects), len(regions), output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, required=True, help="BIDS-like proc folder")
    parser.add_argument("-a", "--suffix", required=True,
                        help="CSV filename suffix, case-insensitive, e.g. parental.csv "
                             "or affectedRegions.csv for the detailed atlas")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Destination CSV")
    parser.add_argument("-g", "--group", type=Path,
                        help="Optional CSV: group-name headers with subject IDs below; "
                             "omit to exclude the group column")
    parser.add_argument("-t", "--translationTXT", dest="translation_txt", type=Path,
                        help="Optional TXT with base label IDs and unprefixed acronyms; "
                             "expand to L_ and R_ (ID + 2000) in memory and translate Label_id")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="After export, interactively select regions and lesion-size units "
                             "and save boxplots as PNG/PDF beside the CSV")
    args = parser.parse_args(argv)
    try:
        summary = export_cohort(args.input, args.suffix, args.output, args.group, args.translation_txt)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    print(f"Saved {summary.sessions} sessions from {summary.subjects} subjects, "
          f"{summary.regions} regions to {summary.output}")
    if args.plot:
        try:
            from affected_regions_plots import plot_csv

            plot_csv(summary.output)
        except (EOFError, KeyboardInterrupt):
            parser.exit(2, f"Plot selection cancelled or input unavailable. CSV retained: {summary.output}\n")
        except (ValueError, OSError, ImportError, RuntimeError) as error:
            parser.exit(2, f"Plotting failed: {error}\nCSV retained: {summary.output}\n")


if __name__ == "__main__":
    main()
