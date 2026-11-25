"""Interactive plots of an exported affected-regions cohort CSV.

Matplotlib is loaded only when drawing. Regional log plots exclude zeros;
lesion-size plots include every observation. Whiskers always show min/max.
"""

import colorsys
import csv
import math
import os
from pathlib import Path
import random
import re
import textwrap
import uuid


TOTAL_COLUMNS = {"mm3": "Total stroke volume [mm3]", "%": "Stroke percent of brain [%]"}
UNITS = {"mm3": " [mm3]", "%": " [%]"}


def column_unit(column):
    return next((unit for unit, suffix in UNITS.items() if column.endswith(suffix)), None)


def read_cohort_csv(path):
    """Read the saved CSV, checking finite nonnegative values for plotting."""
    path = Path(path)
    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream, strict=True)
            header = reader.fieldnames or []
            required = ["subject", "session", *TOTAL_COLUMNS.values()]
            if any(column not in header for column in required):
                raise ValueError("Missing subject, session or total stroke columns.")
            if len(header) != len(set(header)):
                raise ValueError("Duplicate column names.")
            columns = [column for column in header
                       if column_unit(column) and column not in TOTAL_COLUMNS.values()]
            if not columns:
                raise ValueError("No regional measurement columns.")
            rows = []
            samples = set()
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"Line {reader.line_num}: incorrect column count.")
                sample = (row["subject"], row["session"])
                if not all(sample) or sample in samples:
                    raise ValueError(f"Line {reader.line_num}: empty or duplicate subject/session.")
                samples.add(sample)
                for column in [*TOTAL_COLUMNS.values(), *columns]:
                    try:
                        value = float(row[column])
                    except ValueError as error:
                        raise ValueError(f"Line {reader.line_num}, {column}: invalid number.") from error
                    if not math.isfinite(value) or value < 0:
                        raise ValueError(f"Line {reader.line_num}, {column}: "
                                         "plot values must be finite and nonnegative.")
                    row[column] = value
                if "group" in header:
                    row["group"] = row["group"].strip() or "Unassigned"
                rows.append(row)
            if not rows:
                raise ValueError("No samples to plot.")
    except (ValueError, csv.Error) as error:
        raise ValueError(f"{path}: {error}") from error
    return columns, rows


def prompt_unit(description, read_input):
    while True:
        value = read_input(f"{description} (mm3 or %): ").strip().lower().replace("³", "3")
        if value in UNITS:
            return value
        print("Please enter mm3 or %.")


def select_regions(columns, read_input):
    print("Available region columns:")
    for column in columns:
        print(f"  {column}")
    while True:
        selection = read_input("Regions (all or comma-separated full column names): ").strip()
        if selection.lower() == "all":
            while True:
                unit = prompt_unit("Regional measurement", read_input)
                selected = [column for column in columns if column_unit(column) == unit]
                if selected:
                    return selected
                print(f"No regional columns with unit {unit}; choose another measurement.")
        # Region names themselves may contain commas. Only a comma following
        # a complete unit suffix separates two selected column names.
        selected = list(dict.fromkeys(part.strip() for part in re.split(
            r"(?<=\[mm3\])\s*,\s*|(?<=\[%\])\s*,\s*", selection)))
        unknown = [column for column in selected if column not in columns]
        if selection and not unknown:
            return selected
        print("Unknown or empty region columns: " + ", ".join(repr(column) for column in unknown)
              + ". Enter full names including [mm3] or [%].")


def natural_key(value):
    return [int(part) if part.isdigit() else part.casefold() for part in re.split(r"(\d+)", value)]


def group_layout(rows):
    sessions = sorted({row["session"] for row in rows}, key=natural_key)
    groups = {row.get("group", "") for row in rows}
    priority = {"tdcs": 0, "sham": 1, "tdcs sham": 1}
    groups = sorted(groups, key=lambda group: (priority.get(group.casefold(), 2), natural_key(group)))
    return sessions, groups


def group_colors(groups):
    palette = ["#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf"]
    colors = {}
    extra = 0
    for group in groups:
        if group.casefold() in ("tdcs", ""):
            colors[group] = "#0000ee"
        elif group.casefold() in ("sham", "tdcs sham"):
            colors[group] = "#ff0000"
        elif group == "Unassigned":
            colors[group] = "#777777"
        else:
            colors[group] = (palette[extra] if extra < len(palette)
                             else colorsys.hsv_to_rgb((extra * 0.618) % 1, .65, .8))
            extra += 1
    return colors


def session_color(group, color, index, count):
    from matplotlib.colors import to_rgb

    fraction = index / max(1, count - 1)
    if group.casefold() in ("tdcs", ""):
        first, last = to_rgb("#5b51ff"), to_rgb("#000080")
    elif group.casefold() in ("sham", "tdcs sham"):
        first, last = to_rgb("#ff0000"), to_rgb("#a00000")
    else:
        last = to_rgb(color)
        first = tuple(.4 + .6 * component for component in last)
    return tuple(start + fraction * (end - start) for start, end in zip(first, last))


def sample_values(rows, column, session, group):
    return [row[column] for row in rows
            if row["session"] == session and row.get("group", "") == group]


def box_statistics(values, positive_only=False):
    from matplotlib.cbook import boxplot_stats

    values = [value for value in values if value > 0] if positive_only else values
    return boxplot_stats(values, whis=(0, 100))[0] if values else None


def pyplot():
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    return plt


def style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_linewidth(2)
    ax.tick_params(axis="both", which="major", width=2, length=7, labelsize=12)
    ax.tick_params(axis="y", which="minor", width=1.2, length=4)


def draw_box(ax, stats, position, width, color):
    return ax.bxp([stats], positions=[position], widths=width, patch_artist=True,
                  showfliers=False, manage_ticks=False,
                  boxprops={"facecolor": color, "edgecolor": "black", "linewidth": 1.2},
                  whiskerprops={"color": "black", "linewidth": 1.2},
                  capprops={"color": "black", "linewidth": 1.2},
                  medianprops={"color": "black", "linewidth": 1.7})


def region_figure(rows, columns, unit):
    plt = pyplot()
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter, NullFormatter

    sessions, groups = group_layout(rows)
    colors = group_colors(groups)
    combinations = [(session, group) for session in sessions for group in groups]
    width = .8 / len(combinations)
    fig, ax = plt.subplots(figsize=(max(8, 1.25 * len(columns)), 6.6))
    excluded = 0
    positive_values = []
    empty_regions = []
    for index, column in enumerate(columns):
        ax.axvspan(index - .49, index + .49, color="#d3d3d3", zorder=0)
        has_positive = False
        for combo_index, (session, group) in enumerate(combinations):
            values = sample_values(rows, column, session, group)
            zeros = sum(value == 0 for value in values)
            excluded += zeros
            print(f"INFO: {column} / {group + ' / ' if group else ''}{session}: "
                  f"{len(values) - zeros} positive values; excluded zeros: {zeros}.")
            stats = box_statistics(values, positive_only=True)
            if stats is None:
                continue
            has_positive = True
            positive_values.extend(value for value in values if value > 0)
            position = index + (combo_index - (len(combinations) - 1) / 2) * width
            color = session_color(group, colors[group], sessions.index(session), len(sessions))
            draw_box(ax, stats, position, width * .86, color)
            if stats["whislo"] == stats["whishi"]:
                # A collapsed box has no visible fill; retain its group/session
                # color even for a single observation or identical observations.
                ax.scatter([position], [stats["med"]], color=[color], edgecolors="black",
                           linewidths=.4, s=18, zorder=4)
        if not has_positive:
            empty_regions.append(index)
    if not positive_values:
        plt.close(fig)
        print(f"WARNING: Skipping {unit} regional plot: no positive values "
              f"(excluded zeros: {excluded}).")
        return None
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    lower = 10 ** math.floor(math.log10(min(positive_values)) - .05)
    upper = 10 ** math.ceil(math.log10(max(positive_values)) + .05)
    ax.set_ylim(lower, upper)
    for index in empty_regions:
        ax.text(index, .04, "No positive\nvalues", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9)
    names = [column[:-len(UNITS[unit])] for column in columns]
    long_names = any(len(name) > 20 for name in names)
    labels = ["\n".join(textwrap.wrap(name, 26)) for name in names]
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(labels, rotation=40 if long_names else 0,
                       ha="right" if long_names else "center", fontweight="bold")
    ax.set_xlim(-.55, len(columns) - .45)
    ax.set_ylabel("% infarct area" if unit == "%" else "mm³ infarct area",
                  fontsize=17, fontweight="bold")
    handles = [Patch(facecolor=session_color(group, colors[group], sessions.index(session), len(sessions)),
                     edgecolor="black", label=" ".join(filter(None, (group, session.removeprefix("ses-")))))
               for session, group in combinations]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1, 1.01),
              ncol=min(4, len(handles)), frameon=False, fontsize=12)
    style_axes(ax)
    fig.subplots_adjust(left=min(.13, 1.05 / fig.get_figwidth()), right=.99,
                        bottom=.30 if long_names else .17, top=.80)
    fig.text(.5, .025, f"Positive values only; excluded zeros: {excluded}. "
             "Boxes: Q1–Q3; median; whiskers: min–max.", ha="center", fontsize=10)
    return fig


def lesion_figure(rows, unit):
    plt = pyplot()
    from matplotlib.patches import Patch

    sessions, groups = group_layout(rows)
    colors = group_colors(groups)
    width = .7 / len(groups)
    fig, ax = plt.subplots(figsize=(max(6, 1.8 * len(sessions)), 6))
    jitter = random.Random(0)
    for index, session in enumerate(sessions):
        for group_index, group in enumerate(groups):
            values = sample_values(rows, TOTAL_COLUMNS[unit], session, group)
            stats = box_statistics(values)
            if stats is None:
                continue
            position = index + (group_index - (len(groups) - 1) / 2) * width
            draw_box(ax, stats, position, width * .7, colors[group])
            points = [position + jitter.uniform(-width * .22, width * .22) for _ in values]
            ax.scatter(points, values, color=colors[group], edgecolors="black",
                       linewidths=.3, s=35, zorder=4)
    ax.set_xticks(range(len(sessions)))
    ax.set_xticklabels([session.removeprefix("ses-") for session in sessions], fontweight="bold")
    ax.set_xlim(-.6, len(sessions) - .4)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Time point", fontsize=17, fontweight="bold", labelpad=12)
    ax.set_ylabel("Stroke volume (mm³)" if unit == "mm3" else "Stroke percent of brain (%)",
                  fontsize=17, fontweight="bold")
    if groups != [""]:
        ax.legend(handles=[Patch(facecolor=colors[group], edgecolor="black", label=group) for group in groups],
                  loc="lower right", bbox_to_anchor=(1, 1.01), ncol=min(4, len(groups)),
                  frameon=False, fontsize=13)
    style_axes(ax)
    fig.subplots_adjust(left=.17, right=.98, bottom=.18, top=.81)
    fig.text(.5, .025, "All values included. Boxes: Q1–Q3; median; whiskers: min–max.",
             ha="center", fontsize=9)
    return fig


def save_figure(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    width, height = fig.get_size_inches()
    dpi = min(300, 28000 / max(width, height), math.sqrt(80_000_000 / (width * height)))
    if dpi < 300:
        print(f"INFO: Wide PNG rendered at {dpi:.1f} dpi to limit raster size; PDF remains vector.")
    paths = []
    for extension in ("png", "pdf"):
        path = directory / f"{name}.{extension}"
        temporary = directory / f".{name}_{uuid.uuid4().hex}.{extension}"
        try:
            fig.savefig(temporary, dpi=dpi, bbox_inches="tight", facecolor="white")
            os.replace(temporary, path)
        finally:
            if temporary.exists():
                temporary.unlink()
        paths.append(path)
        print(f"Saved plot: {path}")
    return paths


def plot_csv(csv_path, read_input=None):
    """Prompt for selections, then save PNG/PDF plots beside the CSV."""
    columns, rows = read_cohort_csv(csv_path)
    read_input = input if read_input is None else read_input
    selected = select_regions(columns, read_input)
    lesion_unit = prompt_unit("Lesion-size measurement", read_input)
    print("INFO: Box = 25th–75th percentile; line = median; whiskers = minimum to maximum.")
    print("INFO: Regional log plots use positive values only. Lesion-size plots include zeros.")
    print("Selected region columns: " + ", ".join(selected))
    print(f"Selected lesion-size measurement: {lesion_unit}")
    try:
        plt = pyplot()
    except ImportError as error:
        raise ImportError("Matplotlib is required for --plot.") from error
    csv_path = Path(csv_path)
    directory = csv_path.parent / f"{csv_path.stem}_plots"
    paths = []
    with plt.rc_context({"font.family": "DejaVu Sans", "text.usetex": False}):
        for unit in UNITS:
            subset = [column for column in selected if column_unit(column) == unit]
            if not subset:
                continue
            fig = region_figure(rows, subset, unit)
            if fig is not None:
                try:
                    paths.extend(save_figure(fig, directory, f"regions_{'percent' if unit == '%' else 'mm3'}"))
                finally:
                    plt.close(fig)
        fig = lesion_figure(rows, lesion_unit)
        try:
            paths.extend(save_figure(fig, directory, f"lesion_size_{'percent' if lesion_unit == '%' else 'mm3'}"))
        finally:
            plt.close(fig)
    return paths
