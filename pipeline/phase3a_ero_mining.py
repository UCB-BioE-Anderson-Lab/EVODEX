"""
Phase 3a: EVODEX-E operator mining.

This phase takes the filtered EVODEX-P reactions (reaction-level entries with
formula differences) and extracts EVODEX-E operators:

- Each EVODEX-E operator is an abstraction of the reaction mechanism around
  the reactive core (using the EVODEX "E" abstraction level).
- Operators are deduplicated by reaction_hash but all supporting EVODEX-P
  sources are retained.

Outputs:
- data/processed/evodex_e_phase3a_all.csv
- data/errors/phase3a_ero_errors.csv
- data/errors/phase3a_evodex_report.txt
"""

import csv
import os
import sys
import time
from collections import defaultdict
from typing import Dict, List

from evodex.operators import extract_operator_by_abstraction
from evodex.utils import reaction_hash
from pipeline.version import __version__

# Allow very large CSV fields (needed for long SMIRKS strings)
csv.field_size_limit(sys.maxsize)


# ---------------------------------------------------------------------------
# Path configuration (simple hard coded relative paths)
# ---------------------------------------------------------------------------

# Project root (.. from pipeline/)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
ERRORS_DIR = os.path.join(DATA_DIR, "errors")

# Phase 2 product (input to this script)
EVODEX_P_FILTERED = os.path.join(PROCESSED_DIR, "evodex_p_filtered.csv")

# Phase 3 outputs
EVODEX_E_PHASE3_ALL = os.path.join(PROCESSED_DIR, "evodex_e_phase3a_all.csv")

# Error log and report
PHASE3_ERROR_LOG = os.path.join(ERRORS_DIR, "phase3a_ero_errors.csv")
PHASE3_REPORT = os.path.join(ERRORS_DIR, "phase3a_evodex_report.txt")


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def ensure_directories() -> None:
    """Ensure that all necessary directories exist."""
    for path in [PROCESSED_DIR, ERRORS_DIR]:
        os.makedirs(path, exist_ok=True)


# ---------------------------------------------------------------------------
# Main operator extraction
# ---------------------------------------------------------------------------


def extract_evodex_e_operators(
    p_filtered_file: str,
    e_output_file: str,
    error_log_file: str,
) -> Dict[str, Dict]:
    """
    Extract EVODEX-E operators from EVODEX-P filtered reactions.

    Parameters
    ----------
    p_filtered_file:
        Path to EVODEX-P filtered CSV (id, smirks, formula_diff, sources).
    e_output_file:
        Path to write EVODEX-E operators (id, smirks, sources).
    error_log_file:
        Path to write extraction error log.

    Returns
    -------
    operator_map:
        Mapping from operator hash to a dict with keys:
        - "smirks": operator SMIRKS
        - "sources": set of EVODEX-P ids supporting this operator
    """
    # Load EVODEX-P filtered entries
    p_entries: List[Dict] = []

    with open(p_filtered_file, "r") as pfile:
        reader = csv.DictReader(pfile)
        for row in reader:
            p_entries.append(row)

    # Extract EVODEX-E operators and record supporting sources
    operator_map: Dict[str, Dict] = defaultdict(lambda: {"smirks": None, "sources": set()})

    with open(error_log_file, "w", newline="") as errorfile:
        err_writer = csv.DictWriter(
            errorfile, fieldnames=["id", "smirks", "error_message"]
        )
        err_writer.writeheader()

        for i, row in enumerate(p_entries, start=1):
            if i % 100 == 0:
                print(f"Processing EVODEX-P row {i} for operator extraction...")

            rxn_id = row.get("id", "")
            smirks = row.get("smirks", "")

            if not smirks:
                err_writer.writerow(
                    {
                        "id": rxn_id,
                        "smirks": smirks,
                        "error_message": "Missing reaction SMIRKS.",
                    }
                )
                continue

            try:
                # EVODEX-E operator: use the "E" abstraction level
                op_smirks = extract_operator_by_abstraction(smirks, "E")

                # Skip degenerate operators that are empty on one side
                if op_smirks.startswith(">>") or op_smirks.endswith(">>"):
                    continue

                op_hash = reaction_hash(op_smirks)
                operator_map[op_hash]["smirks"] = op_smirks
                operator_map[op_hash]["sources"].add(rxn_id)

            except Exception as exc:
                err_writer.writerow(
                    {
                        "id": rxn_id,
                        "smirks": smirks,
                        "error_message": str(exc),
                    }
                )

    # Write all extracted EVODEX-E operators to output file (all sources retained)
    with open(e_output_file, "w", newline="") as outfile:
        fieldnames = ["id", "smirks", "sources"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for op_hash, data in operator_map.items():
            selected_sources = sorted(data["sources"]) if data["sources"] else []
            writer.writerow(
                {
                    "id": op_hash,
                    "smirks": data["smirks"],
                    "sources": ",".join(selected_sources),
                }
            )

    return operator_map


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
    start_time = time.time()
    print("Phase 3 EVODEX-E mining started...")

    ensure_directories()

    print(f"Reading EVODEX-P filtered from: {EVODEX_P_FILTERED}")
    print(f"Writing EVODEX-E operators to: {EVODEX_E_PHASE3_ALL}")

    operator_map = extract_evodex_e_operators(
        EVODEX_P_FILTERED,
        EVODEX_E_PHASE3_ALL,
        PHASE3_ERROR_LOG,
    )

    # Statistics on extracted operators
    import statistics

    num_sources_list = [len(data["sources"]) for data in operator_map.values()]
    min_sources = min(num_sources_list) if num_sources_list else 0
    max_sources = max(num_sources_list) if num_sources_list else 0
    mean_sources = statistics.mean(num_sources_list) if num_sources_list else 0
    median_sources = statistics.median(num_sources_list) if num_sources_list else 0

    print("Statistics for extracted EVODEX-E operators:")
    print(f"  Number of operators: {len(operator_map)}")
    print(
        f"  Sources per operator: min={min_sources}, "
        f"max={max_sources}, mean={mean_sources:.2f}, median={median_sources}"
    )

    # Save report
    with open(PHASE3_REPORT, "w") as report_file:
        report_file.write(f"EVODEX-E Operator Extraction Report - Phase 3 Mining (v{__version__})\n")
        report_file.write(f"Number of operators: {len(operator_map)}\n")
        report_file.write("Sources per operator:\n")
        report_file.write(f"  min: {min_sources}\n")
        report_file.write(f"  max: {max_sources}\n")
        report_file.write(f"  mean: {mean_sources:.2f}\n")
        report_file.write(f"  median: {median_sources}\n")

    # Report extraction errors
    try:
        with open(PHASE3_ERROR_LOG, "r") as errfile:
            # subtract header
            error_count = max(sum(1 for _ in errfile) - 1, 0)
    except Exception:
        error_count = 0

    if error_count == 0:
        print("✓ No errors encountered during operator extraction.")
    else:
        print(f"Number of extraction errors recorded: {error_count}")

    elapsed_time = time.time() - start_time
    print(f"Phase 3 EVODEX-E mining completed in {elapsed_time:.2f} seconds.")


if __name__ == "__main__":
    main()