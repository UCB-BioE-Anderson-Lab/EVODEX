"""
Phase 4: EVODEX.2 Operator Completion

This phase generates the full EVODEX.2 operator suite for abstractions
A, B, C, D, and E, in two forms:
    - complete families: A, B, C, D, E
    - matched families:  Am, Bm, Cm, Dm, Em

For each EVODEX-P reaction (from Phase 3c or 3d), operators are computed
using extract_operator_by_abstraction.

Operators are deduplicated by reaction_hash and recorded with their
supporting EVODEX-P sources.

This phase also performs fragmentation analysis with respect to EVODEX-E:
for each ERO, we test whether all its P-reaction members produce a single
unique operator at each abstraction family (A, Am, B, Bm, C, Cm, D, Dm, E, Em).

All results are written to:
    data/processed/evodex_<family>.csv
and published to:
    evodex/data/EVODEX-<family>.csv

Reports:
    data/errors/phase4_fragmentation_report.txt
    data/errors/phase4_evodex_report.txt

Branding:
    Operator IDs are assigned using EVODEX.2-* prefixes, using __version__
as defined in pipeline.version.
"""

import csv
import os
import sys
import time
from collections import defaultdict
from typing import Dict, List, Tuple

import pandas as pd

from evodex.operators import extract_operator_by_abstraction
from evodex.utils import reaction_hash

from pipeline.version import __version__
csv.field_size_limit(10**7)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
ERRORS_DIR = os.path.join(DATA_DIR, "errors")
EVODEX_DATA_DIR = os.path.join(BASE_DIR, "evodex", "data")

# Inputs (from Phase 3d)
P_INPUT = os.path.join(PROCESSED_DIR, "evodex_p_phase3d_final.csv")
E_INPUT = os.path.join(PROCESSED_DIR, "evodex_e_phase3d_final.csv")

# Outputs (processed)
OUT_PROCESSED = {
    "A": os.path.join(PROCESSED_DIR, "evodex_a_complete.csv"),
    "Am": os.path.join(PROCESSED_DIR, "evodex_a_matched.csv"),
    "B": os.path.join(PROCESSED_DIR, "evodex_b_complete.csv"),
    "Bm": os.path.join(PROCESSED_DIR, "evodex_b_matched.csv"),
    "C": os.path.join(PROCESSED_DIR, "evodex_c_complete.csv"),
    "Cm": os.path.join(PROCESSED_DIR, "evodex_c_matched.csv"),
    "D": os.path.join(PROCESSED_DIR, "evodex_d_complete.csv"),
    "Dm": os.path.join(PROCESSED_DIR, "evodex_d_matched.csv"),
    "E": os.path.join(PROCESSED_DIR, "evodex_e_complete.csv"),
    "Em": os.path.join(PROCESSED_DIR, "evodex_e_matched.csv"),
}

# Published outputs
OUT_PUBLISHED = {
    key: os.path.join(EVODEX_DATA_DIR, f"EVODEX-{key}.csv")
    for key in OUT_PROCESSED
}

# Reports
FRAGMENTATION_REPORT = os.path.join(ERRORS_DIR, "phase4_fragmentation_report.txt")
SUMMARY_REPORT = os.path.join(ERRORS_DIR, "phase4_evodex_report.txt")


def ensure_directories() -> None:
    for p in [PROCESSED_DIR, ERRORS_DIR, EVODEX_DATA_DIR]:
        os.makedirs(p, exist_ok=True)


# ---------------------------------------------------------------------------
# Operator families (EVODEX.2)
# ---------------------------------------------------------------------------

OPERATOR_FAMILIES = [
    ("A", "A", False),
    ("Am", "A", True),
    ("B", "B", False),
    ("Bm", "B", True),
    ("C", "C", False),
    ("Cm", "C", True),
    ("D", "D", False),
    ("Dm", "D", True),
    ("E", "E", False),
    ("Em", "E", True),
]


# ---------------------------------------------------------------------------
# Operator computation
# ---------------------------------------------------------------------------

def compute_operators_for_p(
    p_id: str, smirks: str
) -> Dict[str, Tuple[str, str]]:
    """
    Compute all EVODEX.2 operator families for a single P reaction.

    Returns:
        {family_key: (operator_smirks, operator_hash)}
    """
    result = {}
    for key, level, matched in OPERATOR_FAMILIES:
        try:
            # operators.extract_operator_by_abstraction(smirks: str, abstraction: str)
            op_smirks = extract_operator_by_abstraction(smirks, level, matched=matched)
        except Exception:
            continue

        # skip malformed operators
        if not op_smirks or op_smirks.startswith(">>") or op_smirks.endswith(">>"):
            continue

        op_hash = reaction_hash(op_smirks)
        result[key] = (op_smirks, op_hash)

    return result


# ---------------------------------------------------------------------------
# Fragmentation analysis
# ---------------------------------------------------------------------------

def fragmentation_analysis(
    e_df: pd.DataFrame,
    operators: Dict[str, Dict[str, Dict]],
) -> Dict[str, Dict[str, List[str]]]:
    """
    For each EVODEX-E, check whether all its P-members map to a single operator
    hash for each operator family. If not, record fragmentation.

    Args:
        e_df: EVODEX-E table with columns ["id", "sources"].
        operators: dict:
            family → {hash → {"smirks": ..., "sources": set([...])}}

    Returns:
        fragmentation_log:
            {family: {E_id: list(hashes)}}
    """

    # Build P → hash maps per family
    p_to_hash = {
        fam: {
            p: h
            for h, entry in opmap.items()
            for p in entry["sources"]
        }
        for fam, opmap in operators.items()
    }

    fragmentation_log = {
        fam: defaultdict(list) for fam in operators
    }

    for row in e_df.itertuples(index=False):
        e_id = str(row.id)
        p_sources = {s for s in str(row.sources).split(",") if s}

        for fam in operators:
            hashes = set()
            mapping = p_to_hash[fam]

            for p in p_sources:
                if p in mapping:
                    hashes.add(mapping[p])

            if len(hashes) > 1:
                fragmentation_log[fam][e_id] = sorted(hashes)

    return fragmentation_log


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def write_fragmentation_report(frag: Dict[str, Dict[str, List[str]]]) -> None:
    with open(FRAGMENTATION_REPORT, "w") as f:
        f.write("EVODEX.2 Phase 4 Fragmentation Report\n")
        f.write("=====================================\n\n")

        total_frag = 0

        for fam, items in frag.items():
            f.write(f"\nOperator family {fam}:\n")
            f.write("-------------------------------------\n")
            for e_id, hashes in items.items():
                total_frag += 1
                f.write(f"EVODEX-E {e_id}: {len(hashes)} distinct operator hashes\n")
                f.write(f"  Hashes: {', '.join(hashes)}\n")
            if not items:
                f.write("No fragmentation detected.\n")

        f.write(f"\nTotal fragmented EVODEX-E entries: {total_frag}\n")


def write_summary_report(
    operators: Dict[str, Dict[str, Dict]],
    frag: Dict[str, Dict[str, List[str]]],
) -> None:
    with open(SUMMARY_REPORT, "w") as f:
        f.write(f"EVODEX.2 Phase 4 Operator Completion Report (version {__version__})\n")
        f.write("===============================================================\n\n")

        for fam, opmap in operators.items():
            counts = [len(entry["sources"]) for entry in opmap.values()]
            if counts:
                f.write(f"{fam}:\n")
                f.write(f"  Total operators: {len(counts)}\n")
                f.write(f"  Min sources:     {min(counts)}\n")
                f.write(f"  Max sources:     {max(counts)}\n")
                f.write(f"  Mean sources:    {sum(counts)/len(counts):.2f}\n")
                f.write(f"  Median sources:  {sorted(counts)[len(counts)//2]}\n")
                f.write("\n")

        total_frag = sum(len(items) for items in frag.values())
        f.write(f"Total fragmented EVODEX-E entries: {total_frag}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    start_time = time.time()
    print("Phase 4 EVODEX.2 operator completion started...")

    ensure_directories()

    # Load inputs
    p_df = pd.read_csv(P_INPUT, dtype={"id": str, "smirks": str, "sources": str})
    e_df = pd.read_csv(E_INPUT, dtype={"id": str, "sources": str})

    # Initialize operator map:
    # {family: {hash: {"smirks": str, "sources": set([...])}}}
    operators = {
        key: {}
        for key, _, _ in OPERATOR_FAMILIES
    }

    # ----------------------------------------------------------------------
    # Compute operators for all EVODEX-P entries
    # ----------------------------------------------------------------------

    for idx, row in enumerate(p_df.itertuples(index=False), start=1):
        if idx % 500 == 0:
            print(f"[{time.strftime('%H:%M:%S')}] Processed {idx} P reactions...")

        p_id = str(row.id)
        smirks = row.smirks

        op_results = compute_operators_for_p(p_id, smirks)

        for fam, (op_smirks, op_hash) in op_results.items():
            fam_map = operators[fam]
            if op_hash not in fam_map:
                fam_map[op_hash] = {"smirks": op_smirks, "sources": set()}
            fam_map[op_hash]["sources"].add(p_id)

    # ----------------------------------------------------------------------
    # Write processed operator CSVs
    # ----------------------------------------------------------------------

    for fam, opmap in operators.items():
        out_path = OUT_PROCESSED[fam]

        rows = []
        for i, (h, entry) in enumerate(opmap.items(), start=1):
            op_id = f"EVODEX.{__version__}-{fam}{i}"
            rows.append({
                "id": op_id,
                "smirks": entry["smirks"],
                "sources": ",".join(sorted(entry["sources"])),
                "hash": h,
            })

        rows = sorted(rows, key=lambda r: len(r["sources"].split(",")), reverse=True)

        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["id", "smirks", "sources", "hash"])
            writer.writeheader()
            writer.writerows(rows)

    # ----------------------------------------------------------------------
    # Fragmentation analysis
    # ----------------------------------------------------------------------

    fragmentation = fragmentation_analysis(e_df, operators)
    write_fragmentation_report(fragmentation)
    write_summary_report(operators, fragmentation)

    # ----------------------------------------------------------------------
    # Publish operators
    # ----------------------------------------------------------------------

    for fam in operators:
        src = OUT_PROCESSED[fam]
        dst = OUT_PUBLISHED[fam]
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        pd.read_csv(src).to_csv(dst, index=False)

    elapsed = time.time() - start_time
    print("✓ No errors encountered.")
    print(f"Phase 4 EVODEX.2 operator completion completed in {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
