#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import sys
from typing import List, Dict, Any

import matplotlib.pyplot as plt

csv.field_size_limit(sys.maxsize)


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _parse_source_count(sources_field: str) -> int:
    if sources_field is None:
        return 0
    s = sources_field.strip()
    if not s:
        return 0
    return len([p for p in s.split(",") if p.strip()])


def count_sources(input_csv: str) -> List[Dict[str, Any]]:
    entries: List[Dict[str, Any]] = []
    with open(input_csv, "r", newline="") as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            entries.append(
                {
                    "id": row["id"],
                    "source_count": _parse_source_count(row.get("sources", "")),
                }
            )
    entries.sort(key=lambda x: x["source_count"], reverse=True)
    return entries


def write_source_counts(entries: List[Dict[str, Any]], output_csv: str) -> None:
    with open(output_csv, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["id", "source_count"])
        writer.writeheader()
        for e in entries:
            writer.writerow(e)


def plot_ranked_support_svg_data_only(entries: List[Dict[str, Any]], output_svg: str) -> None:
    counts = [e["source_count"] for e in entries if e["source_count"] > 0]
    if not counts:
        raise RuntimeError("No nonzero source counts found; cannot plot.")

    n = len(counts)
    x = list(range(n))
    y = counts

    # No text in the figure. Avoid any font issues by not drawing any text.
    fig = plt.figure(figsize=(5.2, 2.6), facecolor="white")

    # Axes fills the whole canvas
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor("white")

    # Log scale (data-only; no ticks/labels)
    ax.set_yscale("log")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(1, max(y))

    # Filled step polygon (single object) from y down to baseline=1
    ax.fill_between(
        x,
        y,
        1,
        step="pre",
        facecolor="black",
        edgecolor="none",
        linewidth=0,
    )

    # Hide everything except the frame
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # Keep a rectangular border (black)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)

    # Absolutely no grid
    ax.grid(False)

    # Save tightly so the SVG is just the panel
    fig.savefig(output_svg, format="svg", bbox_inches="tight", pad_inches=0, facecolor="white")
    plt.close(fig)


def main() -> None:
    base_dir = os.path.join("analysis", "unfiltered_mining", "output")
    _ensure_dir(base_dir)

    input_csv = os.path.join(base_dir, "unfiltered_eros.csv")
    output_counts_csv = os.path.join(base_dir, "ero_source_counts.csv")
    output_svg = os.path.join(base_dir, "unfiltered_ero_support_rank.svg")

    entries = count_sources(input_csv)
    write_source_counts(entries, output_counts_csv)
    plot_ranked_support_svg_data_only(entries, output_svg)

    print(f"Wrote: {output_counts_csv}")
    print(f"Wrote: {output_svg}")
    print(f"Operators counted: {len(entries)}")


if __name__ == "__main__":
    main()
