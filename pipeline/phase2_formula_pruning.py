"""
Phase 2: EVODEX formula pruning.

This script takes EVODEX-R preliminary reactions and derives:

- EVODEX-P (preliminary): reaction level, annotated with overall formula
  change between reactants and products.
- EVODEX-F (preliminary): formula level, grouping reactions by the same
  net formula change.
- EVODEX-F (filtered): formula entries with support above a threshold.
- EVODEX-P (filtered): reaction entries whose formula group passed the
  support threshold.

Schema notes (current):
- EVODEX-F:
    id           = canonical formula key string (e.g., "C+2;H-3;O+1" or "0")
    formula_diff = JSON-serialized dict (e.g., {"C": 2, "H": -3, "O": 1})
- EVODEX-P filtered:
    formula_diff = formula key string (matches EVODEX-F id)

Outputs:

- data/processed/temp_evodex_p_preliminary.csv
- data/processed/temp_evodex_f_preliminary.csv
- data/processed/evodex_f_filtered.csv
- data/processed/evodex_p_filtered.csv
- data/errors/split_reactions_errors.csv
- data/errors/formula_errors.csv
- data/errors/Phase2_evodex_report.txt
"""

import csv
import json
import os
import time
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Tuple, Set

# Allow very large SMIRKS / sources fields
csv.field_size_limit(10**7)

from rdkit import Chem
from rdkit.Chem import AllChem

from evodex.splitting import split_reaction
from evodex.utils import reaction_hash
from pipeline.version import __version__


# ---------------------------------------------------------------------------
# Path configuration (simple hard coded relative paths)
# ---------------------------------------------------------------------------

# Project root (.. from pipeline/)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
ERRORS_DIR = os.path.join(DATA_DIR, "errors")

# Phase 1 product (input to this script)
EVODEX_R_PRELIMINARY = os.path.join(PROCESSED_DIR, "temp_evodex_r_preliminary.csv")

# Phase 2 outputs
EVODEX_P_PRELIMINARY = os.path.join(PROCESSED_DIR, "temp_evodex_p_preliminary.csv")
EVODEX_F_PRELIMINARY = os.path.join(PROCESSED_DIR, "temp_evodex_f_preliminary.csv")
EVODEX_F_FILTERED = os.path.join(PROCESSED_DIR, "evodex_f_filtered.csv")
EVODEX_P_FILTERED = os.path.join(PROCESSED_DIR, "evodex_p_filtered.csv")

# Error logs and report
SPLIT_REACTIONS_ERRORS = os.path.join(ERRORS_DIR, "split_reactions_errors.csv")
FORMULA_ERRORS = os.path.join(ERRORS_DIR, "formula_errors.csv")
PHASE2_REPORT = os.path.join(ERRORS_DIR, "Phase2_evodex_report.txt")

# Minimum number of reactions required to keep a formula group
MIN_FORMULA_SUPPORT = 5


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def ensure_directories() -> None:
    for path in [PROCESSED_DIR, ERRORS_DIR]:
        os.makedirs(path, exist_ok=True)


def log_progress(prefix: str, idx: int, step: int = 100) -> None:
    if idx % step == 0:
        print(f"[{time.strftime('%H:%M:%S')}] {prefix} {idx}...")


# ---------------------------------------------------------------------------
# Chemistry helpers
# ---------------------------------------------------------------------------


def split_reaction_smiles(rxn_smiles: str) -> Tuple[List[Chem.Mol], List[Chem.Mol]]:
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
    except Exception as exc:
        raise ValueError(f"Could not parse reaction SMILES: {rxn_smiles}") from exc

    reactant_mols = list(rxn.GetReactants())
    product_mols = list(rxn.GetProducts())

    if not reactant_mols and not product_mols:
        raise ValueError(f"Reaction has no reactants and no products: {rxn_smiles}")

    return reactant_mols, product_mols


def _composition_from_mols(mols: Iterable[Chem.Mol]) -> Counter:
    comp: Counter = Counter()
    for mol in mols:
        for atom in mol.GetAtoms():
            comp[atom.GetSymbol()] += 1
    return comp


def _format_formula_diff(diff: Dict[str, int]) -> str:
    if not diff:
        return "0"

    parts: List[str] = []
    for elem in sorted(diff):
        delta = int(diff[elem])
        if delta == 0:
            continue
        sign = "+" if delta > 0 else ""
        parts.append(f"{elem}{sign}{delta}")
    return ";".join(parts) if parts else "0"


def calculate_formula_diff_obj(rxn_smiles: str) -> Dict[str, int]:
    """
    Compute the net formula change between reactants and products.

    Returns a dict like {"C": 2, "H": -3, "O": 1}.
    """
    reactants, products = split_reaction_smiles(rxn_smiles)
    comp_r = _composition_from_mols(reactants)
    comp_p = _composition_from_mols(products)

    all_elems = set(comp_r) | set(comp_p)
    diff: Dict[str, int] = {}
    for elem in all_elems:
        delta = int(comp_p[elem] - comp_r[elem])
        if delta != 0:
            diff[elem] = delta

    return diff


# ---------------------------------------------------------------------------
# Stage 1: derive EVODEX-P from EVODEX-R
# ---------------------------------------------------------------------------


def process_split_reactions(
    input_file: str,
    p_output_file: str,
    split_error_file: str,
) -> int:
    error_rows: List[Dict] = []
    hash_to_ids: Dict[str, Set[str]] = defaultdict(set)
    hash_to_example_smirks: Dict[str, str] = {}

    with open(input_file, "r") as infile:
        reader = csv.DictReader(infile)
        for idx, row in enumerate(reader, start=1):
            log_progress("split_reactions processing row", idx)

            smirks = row.get("smirks", "")
            rxn_id = row.get("id", "")
            sources = row.get("sources", "")

            if not smirks:
                error_rows.append(
                    {
                        "id": rxn_id,
                        "smirks": smirks,
                        "sources": sources,
                        "error": "Missing reaction SMIRKS.",
                    }
                )
                continue

            try:
                split_reactions = split_reaction(smirks)
            except Exception as exc:
                error_rows.append(
                    {
                        "id": rxn_id,
                        "smirks": smirks,
                        "sources": sources,
                        "error": str(exc),
                    }
                )
                continue

            for reaction in split_reactions:
                try:
                    rxn_hash = reaction_hash(reaction)
                except Exception as exc:
                    error_rows.append(
                        {
                            "id": rxn_id,
                            "smirks": reaction,
                            "sources": sources,
                            "error": f"Failed to hash split reaction: {exc}",
                        }
                    )
                    continue

                hash_to_ids[rxn_hash].add(rxn_id)
                if rxn_hash not in hash_to_example_smirks:
                    hash_to_example_smirks[rxn_hash] = reaction

    p_rows_written = 0
    with open(p_output_file, "w", newline="") as outfile:
        fieldnames = ["id", "smirks", "sources"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for rxn_hash, id_set in hash_to_ids.items():
            example_smirks = hash_to_example_smirks[rxn_hash]
            sources_str = ",".join(sorted(id_set))
            writer.writerow({"id": rxn_hash, "smirks": example_smirks, "sources": sources_str})
            p_rows_written += 1

    if error_rows:
        fieldnames = ["id", "smirks", "sources", "error"]
        with open(split_error_file, "w", newline="") as ef:
            writer = csv.DictWriter(ef, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(error_rows)

    return p_rows_written


# ---------------------------------------------------------------------------
# Stage 2: formula level grouping and pruning
# ---------------------------------------------------------------------------


def process_formula_pruning(
    p_input_file: str,
    f_output_file: str,
    f_filtered_file: str,
    p_filtered_file: str,
    formula_error_file: str,
    min_support: int = MIN_FORMULA_SUPPORT,
) -> Tuple[int, int, int, int]:
    error_rows: List[Dict] = []
    formula_groups: Dict[str, Dict] = {}   # formula_key -> {"formula_json": str, "sources": set}
    formula_support: Counter = Counter()   # formula_key -> count
    p_id_to_formula_key: Dict[str, str] = {}

    # First pass: read EVODEX-P and accumulate groups
    with open(p_input_file, "r") as infile:
        reader = csv.DictReader(infile)
        n_p_rows = 0
        for idx, row in enumerate(reader, start=1):
            log_progress("formula_pruning processing row", idx)
            n_p_rows += 1

            p_id = row.get("id", "")
            smirks = row.get("smirks", "")
            sources = row.get("sources", "")

            if not smirks:
                error_rows.append(
                    {"id": p_id, "smirks": smirks, "formula_diff": "", "sources": sources, "error": "Missing reaction SMIRKS."}
                )
                continue

            try:
                diff_obj = calculate_formula_diff_obj(smirks)       # dict
                formula_key = _format_formula_diff(diff_obj)        # canonical string for grouping + id
                formula_json = json.dumps(diff_obj, sort_keys=True) # serialized dict for storage
            except Exception as exc:
                error_rows.append(
                    {"id": p_id, "smirks": smirks, "formula_diff": "", "sources": sources, "error": str(exc)}
                )
                continue

            if p_id:
                p_id_to_formula_key[p_id] = formula_key

            if formula_key not in formula_groups:
                formula_groups[formula_key] = {"formula_json": formula_json, "sources": set()}

            formula_support[formula_key] += 1
            if p_id:
                formula_groups[formula_key]["sources"].add(p_id)

    # Write EVODEX-F preliminary
    n_f_rows = 0
    with open(f_output_file, "w", newline="") as outfile:
        fieldnames = ["id", "formula_diff", "support", "sources"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for formula_key in sorted(formula_groups.keys()):
            support = int(formula_support[formula_key])
            sources_set = formula_groups[formula_key]["sources"]
            sources_str = ",".join(sorted(sources_set)) if sources_set else ""

            writer.writerow(
                {
                    "id": formula_key,
                    "formula_diff": formula_groups[formula_key]["formula_json"],
                    "support": support,
                    "sources": sources_str,
                }
            )
            n_f_rows += 1

    kept_formula_keys = {k for k, count in formula_support.items() if count >= min_support}

    # Write EVODEX-F filtered
    n_f_kept = 0
    with open(f_filtered_file, "w", newline="") as outfile:
        fieldnames = ["id", "formula_diff", "support", "sources"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        with open(f_output_file, "r") as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                if row["id"] in kept_formula_keys:
                    writer.writerow(row)
                    n_f_kept += 1

    # Write EVODEX-P filtered (only reactions whose formula group was kept)
    n_p_kept = 0
    with open(p_input_file, "r") as infile, open(p_filtered_file, "w", newline="") as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ["id", "smirks", "formula_diff", "sources"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            p_id = row.get("id", "")
            formula_key = p_id_to_formula_key.get(p_id)
            if not formula_key or formula_key not in kept_formula_keys:
                continue

            writer.writerow(
                {
                    "id": p_id,
                    "smirks": row.get("smirks", ""),
                    "formula_diff": formula_key,  # store the canonical key string (matches EVODEX-F id)
                    "sources": row.get("sources", ""),
                }
            )
            n_p_kept += 1

    if error_rows:
        fieldnames = ["id", "smirks", "formula_diff", "sources", "error"]
        with open(formula_error_file, "w", newline="") as ef:
            writer = csv.DictWriter(ef, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(error_rows)

    return n_p_rows, n_f_rows, n_f_kept, n_p_kept


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
    start_time = time.time()
    print("Phase 2 formula pruning started...")

    ensure_directories()

    print(f"Reading EVODEX-R preliminary from: {EVODEX_R_PRELIMINARY}")
    print(f"Writing EVODEX-P preliminary to: {EVODEX_P_PRELIMINARY}")
    p_rows_written = process_split_reactions(
        EVODEX_R_PRELIMINARY,
        EVODEX_P_PRELIMINARY,
        SPLIT_REACTIONS_ERRORS,
    )

    print(f"Processing formula pruning from: {EVODEX_P_PRELIMINARY}")
    n_p_rows, n_f_rows, n_f_kept, n_p_kept = process_formula_pruning(
        EVODEX_P_PRELIMINARY,
        EVODEX_F_PRELIMINARY,
        EVODEX_F_FILTERED,
        EVODEX_P_FILTERED,
        FORMULA_ERRORS,
        MIN_FORMULA_SUPPORT,
    )

    report_lines: List[str] = [
        f"EVODEX Phase 2 Report (version {__version__})",
        "==============================================",
        f"EVODEX-P preliminary reactions read: {n_p_rows}",
        f"EVODEX-P preliminary reactions written: {p_rows_written}",
        f"EVODEX-F preliminary formula groups: {n_f_rows}",
        f"EVODEX-F filtered (support >= {MIN_FORMULA_SUPPORT}): {n_f_kept}",
        f"EVODEX-P filtered reactions: {n_p_kept}",
        "",
        f"Split reactions errors file: {SPLIT_REACTIONS_ERRORS}",
        f"Formula errors file: {FORMULA_ERRORS}",
    ]

    for line in report_lines:
        print(line)

    with open(PHASE2_REPORT, "w") as report_file:
        report_file.write("\n".join(report_lines))

    split_errors_exist = os.path.exists(SPLIT_REACTIONS_ERRORS) and os.path.getsize(SPLIT_REACTIONS_ERRORS) > 0
    formula_errors_exist = os.path.exists(FORMULA_ERRORS) and os.path.getsize(FORMULA_ERRORS) > 0

    if not split_errors_exist and not formula_errors_exist:
        print("✓ No errors encountered.")

    elapsed = time.time() - start_time
    print(f"Phase 2 formula pruning completed in {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()