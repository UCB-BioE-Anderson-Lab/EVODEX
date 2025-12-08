"""
Phase 3b: EVODEX-E pruning.

This phase takes the raw EVODEX-E operators from Phase 3a (evodex_e_phase3_all)
and performs a two-step pruning:

1. Remove operators whose reactants have no mapped atoms (no atom-map numbers).
2. Among the remaining operators, retain only those with at least a minimum
   number of supporting EVODEX-P reactions, and cap the number of sources
   kept per operator.

Outputs
-------
- data/processed/evodex_e_phase3b_preliminary.csv
    Pruned EVODEX-E operators, each with up to MAX_SOURCES_PER_OPERATOR
    supporting EVODEX-P ids.
- data/processed/evodex_p_phase3b_retained.csv
    Subset of EVODEX-P filtered entries whose ids appear in the surviving
    EVODEX-E operators.
- data/errors/phase3b_evodex_report.txt
    Text report with pruning statistics.
"""

import csv
import os
import time
import statistics
from collections import defaultdict
from typing import Dict, List, Set

from rdkit import Chem
from rdkit.Chem import AllChem


from pipeline.version import __version__

# Pruning thresholds
MIN_SOURCES_PER_OPERATOR = 10
MAX_SOURCES_PER_OPERATOR = 5

# Allow very large CSV fields (needed for long SMIRKS strings)
csv.field_size_limit(10**7)

# ---------------------------------------------------------------------------
# Path configuration (simple hard-coded relative paths)
# ---------------------------------------------------------------------------

# Project root (.. from pipeline/)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
ERRORS_DIR = os.path.join(DATA_DIR, "errors")

# Inputs
EVODEX_E_PHASE3_ALL = os.path.join(PROCESSED_DIR, "evodex_e_phase3a_all.csv")
EVODEX_P_FILTERED = os.path.join(PROCESSED_DIR, "evodex_p_filtered.csv")

# Outputs for this phase
EVODEX_E_PHASE3B_PRELIM = os.path.join(
    PROCESSED_DIR, "evodex_e_phase3b_preliminary.csv"
)
EVODEX_P_PHASE3B_RETAINED = os.path.join(
    PROCESSED_DIR, "evodex_p_phase3b_retained.csv"
)

# Report
PHASE3B_REPORT = os.path.join(ERRORS_DIR, "phase3b_evodex_report.txt")


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def ensure_directories() -> None:
    """Ensure that all necessary directories exist."""
    for path in [PROCESSED_DIR, ERRORS_DIR]:
        os.makedirs(path, exist_ok=True)


def log_progress(prefix: str, idx: int, step: int = 100) -> None:
    if idx % step == 0:
        print(f"[{time.strftime('%H:%M:%S')}] {prefix} {idx}...")


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------


def has_mapped_atoms_in_reactants(smarts: str) -> bool:
    """
    Return True if at least one atom in the reactants has a positive atom-map
    number in the given operator SMIRKS/SMARTS string.
    """
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
    except Exception:
        return False

    for reactant in rxn.GetReactants():
        for atom in reactant.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                return True
    return False


def load_evodex_e(filepath: str) -> Dict[str, Dict]:
    """
    Load EVODEX-E operators from CSV.

    Returns a dict keyed by operator id with:
      - 'smirks': operator SMIRKS
      - 'sources': set of EVODEX-P ids
    """
    operators: Dict[str, Dict] = {}

    with open(filepath, "r") as infile:
        reader = csv.DictReader(infile)
        for idx, row in enumerate(reader, start=1):
            log_progress("Reading EVODEX-E", idx)
            op_id = row.get("id", "")
            smirks = row.get("smirks", "")
            sources_str = row.get("sources", "")

            if not op_id or not smirks:
                continue

            sources: Set[str] = set()
            if sources_str:
                # Sources are stored as comma-separated ids
                sources = {s for s in sources_str.split(",") if s}

            operators[op_id] = {
                "smirks": smirks,
                "sources": sources,
            }

    return operators


def prune_evodex_e(operators: Dict[str, Dict]) -> Dict[str, Dict]:
    """
    Apply pruning rules to EVODEX-E operators.

    Steps:
    1. Remove operators with no mapped atoms in reactants.
    2. Keep only operators with at least MIN_SOURCES_PER_OPERATOR sources.
    3. Cap number of sources per operator at MAX_SOURCES_PER_OPERATOR.
    """
    total_operators = len(operators)
    rejected_unmapped = 0
    rejected_low_support = 0

    retained: Dict[str, Dict] = {}

    for idx, (op_id, data) in enumerate(operators.items(), start=1):
        log_progress("Pruning EVODEX-E", idx)
        smirks = data["smirks"]
        sources: Set[str] = data["sources"]

        # Step 1: require at least one mapped atom in the reactants
        if not has_mapped_atoms_in_reactants(smirks):
            rejected_unmapped += 1
            continue

        # Step 2: minimum support threshold
        source_count = len(sources)
        if source_count < MIN_SOURCES_PER_OPERATOR:
            rejected_low_support += 1
            continue

        # Step 3: cap number of sources per operator
        capped_sources = sorted(sources)[:MAX_SOURCES_PER_OPERATOR]

        retained[op_id] = {
            "smirks": smirks,
            "sources": capped_sources,
            "original_source_count": source_count,
        }

    print(f"Total EVODEX-E operators: {total_operators}")
    print(f"Rejected (no mapped atoms in reactants): {rejected_unmapped}")
    print(f"Rejected (support < {MIN_SOURCES_PER_OPERATOR}): {rejected_low_support}")
    print(f"Retained operators: {len(retained)}")

    return retained


def write_pruned_evodex_e(operators: Dict[str, Dict], filepath: str) -> None:
    """
    Write pruned EVODEX-E operators to CSV.
    """
    with open(filepath, "w", newline="") as outfile:
        fieldnames = ["id", "smirks", "sources", "original_source_count"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for op_id, data in operators.items():
            writer.writerow(
                {
                    "id": op_id,
                    "smirks": data["smirks"],
                    "sources": ",".join(data["sources"]),
                    "original_source_count": data["original_source_count"],
                }
            )


def write_retained_evodex_p(
    p_input_file: str,
    p_output_file: str,
    surviving_p_ids: Set[str],
) -> int:
    """
    Write the subset of EVODEX-P filtered rows whose ids are in surviving_p_ids.

    Returns the count of retained P rows.
    """
    retained_count = 0

    with open(p_input_file, "r") as infile, open(
        p_output_file, "w", newline=""
    ) as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames or ["id", "smirks", "formula_diff", "sources"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for idx, row in enumerate(reader, start=1):
            log_progress("Writing retained EVODEX-P", idx)
            rxn_id = row.get("id", "")
            if rxn_id in surviving_p_ids:
                writer.writerow(row)
                retained_count += 1

    return retained_count


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
    start_time = time.time()
    print("Phase 3b EVODEX-E pruning started...")

    ensure_directories()

    print(f"Reading EVODEX-E operators from: {EVODEX_E_PHASE3_ALL}")
    operators = load_evodex_e(EVODEX_E_PHASE3_ALL)

    print("Pruning EVODEX-E operators...")
    pruned_operators = prune_evodex_e(operators)

    print(f"Writing pruned EVODEX-E operators to: {EVODEX_E_PHASE3B_PRELIM}")
    write_pruned_evodex_e(pruned_operators, EVODEX_E_PHASE3B_PRELIM)

    # Collect all P ids that survived in at least one operator
    surviving_p_ids: Set[str] = set()
    for data in pruned_operators.values():
        surviving_p_ids.update(data["sources"])

    print(f"Reading EVODEX-P filtered from: {EVODEX_P_FILTERED}")
    print(f"Writing retained EVODEX-P entries to: {EVODEX_P_PHASE3B_RETAINED}")
    retained_p_count = write_retained_evodex_p(
        EVODEX_P_FILTERED, EVODEX_P_PHASE3B_RETAINED, surviving_p_ids
    )

    # Statistics on sources per operator (using original source counts)
    original_counts = [data["original_source_count"] for data in pruned_operators.values()]
    if original_counts:
        min_sources = min(original_counts)
        max_sources = max(original_counts)
        mean_sources = statistics.mean(original_counts)
        median_sources = statistics.median(original_counts)
    else:
        min_sources = max_sources = mean_sources = median_sources = 0

    print("Statistics for retained EVODEX-E operators (before capping sources):")
    print(f"  Number of operators: {len(pruned_operators)}")
    print(
        f"  Sources per operator: min={min_sources}, "
        f"max={max_sources}, mean={mean_sources:.2f}, median={median_sources}"
    )
    print(f"Retained EVODEX-P reactions: {retained_p_count}")

    # Write report
    report_lines: List[str] = [
        f"EVODEX-E Operator Retention Report (Phase 3b, v{__version__})",
        "===============================================================",
        f"Input EVODEX-E operators: {len(operators)}",
        f"Retained EVODEX-E operators: {len(pruned_operators)}",
        f"Retained EVODEX-P reactions: {retained_p_count}",
        "",
        "Sources per retained EVODEX-E operator (before capping to "
        f"{MAX_SOURCES_PER_OPERATOR}):",
        f"  min: {min_sources}",
        f"  max: {max_sources}",
        f"  mean: {mean_sources:.2f}",
        f"  median: {median_sources}",
    ]

    for line in report_lines:
        print(line)

    with open(PHASE3B_REPORT, "w") as report_file:
        report_file.write("\n".join(report_lines))

    print("✓ No errors encountered.")

    elapsed = time.time() - start_time
    print(f"Phase 3b EVODEX-E pruning completed in {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
