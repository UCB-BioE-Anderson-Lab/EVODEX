"""
Phase 1: EVODEX-R data preparation.

This script prepares the EVODEX-R dataset from raw mapped reactions.
It performs the following steps:

1. Read raw mapped reaction SMILES from data/raw/raw_reactions.csv.
2. Validate that each reaction parses and can be hashed.
3. Write a normalized filtered CSV with basic error logging.
4. Add explicit hydrogens and hydrogen atom maps using hydrogen_mapper.
5. Deduplicate reactions by hash and compute mining statistics.
6. Write:
   - data/processed/hydrogen_mapped_reactions.csv (intermediate)
   - data/processed/temp_evodex_r_preliminary.csv (deduplicated EVODEX-R)
   - data/errors/Phase1_evodex_report.txt (summary report)
   - data/errors/data_preparation_errors.csv (detailed error log, if any).
"""

import csv
import os
import time
import statistics
from collections import Counter
from typing import Dict, List, Tuple

from evodex.hydrogen_mapper import map_hydrogens_in_reaction
from evodex.utils import reaction_hash
from pipeline.version import __version__


# ---------------------------------------------------------------------------
# Path configuration (simple hard-coded relative paths)
# ---------------------------------------------------------------------------

# Project root (.. from pipeline/)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
ERRORS_DIR = os.path.join(DATA_DIR, "errors")

# Phase 1 inputs / outputs
RAW_DATA = os.path.join(RAW_DIR, "raw_reactions.csv")
FILTERED_DATA = os.path.join(PROCESSED_DIR, "filtered_reactions.csv")
HYDROGEN_MAPPED_DATA = os.path.join(PROCESSED_DIR, "hydrogen_mapped_reactions.csv")

# Deduplicated EVODEX-R (Phase 1 product)
EVODEX_R_PRELIMINARY = os.path.join(
    PROCESSED_DIR, "temp_evodex_r_preliminary.csv"
)

# Reports / logs
PHASE1_REPORT = os.path.join(ERRORS_DIR, "Phase1_evodex_report.txt")
ERROR_LOG_CSV = os.path.join(ERRORS_DIR, "data_preparation_errors.csv")


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def ensure_directories() -> None:
    """Ensure that all necessary directories exist."""
    for path in [RAW_DIR, PROCESSED_DIR, ERRORS_DIR]:
        os.makedirs(path, exist_ok=True)


def write_row(writer: csv.DictWriter, data: Dict) -> None:
    writer.writerow(data)


def log_progress(prefix: str, idx: int, step: int = 100) -> None:
    if idx % step == 0:
        print(f"[{time.strftime('%H:%M:%S')}] {prefix} {idx}...")


# ---------------------------------------------------------------------------
# Stage 0: validate and normalize raw reactions
# ---------------------------------------------------------------------------


def process_raw_data(input_file: str, output_file: str) -> Tuple[int, int]:
    """
    Validate raw mapped reaction SMILES and write a normalized CSV.

    Parameters
    ----------
    input_file:
        CSV file containing at least columns "rxn_idx" and "mapped".
    output_file:
        Destination CSV with columns id, smirks, sources, error.

    Returns
    -------
    total_raw, valid
        The number of raw rows processed and the number that passed validation.
    """
    total_raw = 0
    valid = 0

    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ["id", "smirks", "sources", "error"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            # Skip completely empty lines
            if not any(row.values()):
                continue

            total_raw += 1
            log_progress("Processed raw reactions", total_raw)

            rxn_id = row.get("rxn_idx", "")
            mapped = row.get("mapped", "")

            if not mapped:
                new_row = {
                    "id": rxn_id,
                    "smirks": "",
                    "sources": "",
                    "error": "Missing mapped reaction SMILES.",
                }
                write_row(writer, new_row)
                continue

            try:
                # Sanity check that the reaction can be hashed
                reaction_hash(mapped)
                valid += 1
                new_row = {
                    "id": rxn_id,
                    "smirks": mapped,
                    "sources": rxn_id,
                    "error": "",
                }
            except Exception as exc:
                new_row = {
                    "id": rxn_id,
                    "smirks": "",
                    "sources": "",
                    "error": str(exc),
                }

            write_row(writer, new_row)

    return total_raw, valid


# ---------------------------------------------------------------------------
# Generic row-wise processing helper
# ---------------------------------------------------------------------------


def process_data(
    input_file: str,
    output_file: str,
    transformation_function,
    stage_name: str,
    error_log: List[Dict],
) -> None:
    """
    Apply a row-wise transformation function to a CSV file with columns
    id, smirks, sources, error.

    Rows that already have a non-empty error field are passed through
    unchanged. Errors raised by the transformation are captured and
    written to the error log.
    """
    import traceback

    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ["id", "smirks", "sources", "error"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for idx, row in enumerate(reader, start=1):
            log_progress(f"{stage_name} processing row", idx)

            if row.get("error"):
                write_row(writer, row)
                continue

            try:
                new_data = transformation_function(row)
                write_row(writer, new_data)
            except Exception as exc:
                error_msg = f"{exc}\n{traceback.format_exc()}"
                row["error"] = error_msg
                error_log.append({**row, "stage": stage_name})
                write_row(writer, row)


# ---------------------------------------------------------------------------
# Stage 1: add explicit hydrogens and hydrogen atom maps
# ---------------------------------------------------------------------------


def hydrogen_mapping_transform(row: Dict) -> Dict:
    """
    Add explicit hydrogens and hydrogen atom maps to a heavy-atom-mapped
    reaction.

    Parameters
    ----------
    row:
        A dictionary with keys id, smirks, sources, error.

    Returns
    -------
    dict
        A new row with updated smirks and cleared error.
    """
    mapped_with_h = map_hydrogens_in_reaction(row["smirks"])
    return {
        "id": row["id"],
        "smirks": mapped_with_h,
        "sources": row["sources"],
        "error": "",
    }


# ---------------------------------------------------------------------------
# Mining summary and deduplication
# ---------------------------------------------------------------------------


def summarize_mining(mapped_file: str) -> Tuple[List[str], Dict[str, Dict]]:
    """
    Compute EVODEX-R mining statistics and collect unique reactions.

    Parameters
    ----------
    mapped_file:
        CSV file with hydrogen-mapped reactions.

    Returns
    -------
    report_lines:
        Human-readable report lines.
    unique_reactions:
        Mapping from reaction hash to a representative row dict.
    """
    smirks_counts: Counter = Counter()

    with open(mapped_file, "r") as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            if row.get("error"):
                continue
            smirks = row["smirks"]
            rxn_hash = reaction_hash(smirks)
            smirks_counts[rxn_hash] += 1

    unique_ev_r = len(smirks_counts)
    counts = list(smirks_counts.values())

    report_lines: List[str] = [
        f"EVODEX-R Mining Report (version {__version__})",
        "============================================",
        f"Total unique EVODEX-R entries: {unique_ev_r}",
    ]

    if counts:
        report_lines.extend(
            [
                f"Min reactions per EVODEX-R: {min(counts)}",
                f"Max reactions per EVODEX-R: {max(counts)}",
                f"Mean reactions per EVODEX-R: {statistics.mean(counts):.2f}",
                f"Median reactions per EVODEX-R: {statistics.median(counts)}",
            ]
        )
    else:
        report_lines.extend(
            [
                "Min reactions per EVODEX-R: 0",
                "Max reactions per EVODEX-R: 0",
                "Mean reactions per EVODEX-R: 0.00",
                "Median reactions per EVODEX-R: 0",
            ]
        )

    report_lines.append("")

    # Build a representative set of unique reactions
    unique_reactions: Dict[str, Dict] = {}
    with open(mapped_file, "r") as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            if row.get("error"):
                continue
            smirks = row["smirks"]
            rxn_hash = reaction_hash(smirks)
            if rxn_hash not in unique_reactions:
                unique_reactions[rxn_hash] = {
                    "id": rxn_hash,
                    "smirks": smirks,
                    "sources": row["sources"],
                }

    return report_lines, unique_reactions


def write_deduplicated_evodex(unique_reactions: Dict[str, Dict]) -> None:
    """Write the deduplicated EVODEX-R reactions to CSV."""
    with open(EVODEX_R_PRELIMINARY, "w", newline="") as outfile:
        fieldnames = ["id", "smirks", "sources"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for data in unique_reactions.values():
            writer.writerow(data)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
    start_time = time.time()
    print("Phase 1 data preparation started...")

    ensure_directories()

    # Stage 0: validate raw reactions
    print(f"Reading raw reactions from: {RAW_DATA}")
    total_raw, valid_reactions = process_raw_data(RAW_DATA, FILTERED_DATA)

    # Stage 1: hydrogen mapping
    error_log: List[Dict] = []
    print(f"Running hydrogen mapping on: {FILTERED_DATA}")
    process_data(
        FILTERED_DATA,
        HYDROGEN_MAPPED_DATA,
        hydrogen_mapping_transform,
        "hydrogen_mapping",
        error_log,
    )

    # Mining summary and deduplication
    print(f"Summarizing mining statistics from: {HYDROGEN_MAPPED_DATA}")
    report_lines, unique_reactions = summarize_mining(HYDROGEN_MAPPED_DATA)

    percentage_loss = (
        (total_raw - valid_reactions) / total_raw * 100 if total_raw else 0.0
    )

    header_lines = [
        f"Total raw reactions: {total_raw}",
        f"Valid reactions after sanitization: {valid_reactions}",
        f"Compression rate due to reaction deduplication: {percentage_loss:.2f}%",
        (
            f"(i.e., only {100 - percentage_loss:.2f}% of raw reactions are "
            "unique after hashing)"
        ),
        "",
    ]

    full_report = header_lines + report_lines

    for line in full_report:
        print(line)

    # Write report
    with open(PHASE1_REPORT, "w") as report_file:
        report_file.write("\n".join(full_report))

    # Write deduplicated EVODEX-R
    write_deduplicated_evodex(unique_reactions)

    # Write error log if needed
    if error_log:
        fieldnames = sorted(error_log[0].keys())
        with open(ERROR_LOG_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(error_log)
    if not error_log:
        print("✓ No errors encountered.")

    elapsed = time.time() - start_time
    print(f"Phase 1 data preparation completed in {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
