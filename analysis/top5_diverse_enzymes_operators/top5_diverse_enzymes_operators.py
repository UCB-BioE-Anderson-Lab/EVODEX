#!/usr/bin/env python3
"""
Top-5 diverse enzymes: map reactions and extract Bm, Cm, Em operators.

This script:
  1) Chooses 5 enzymes with the largest reaction counts, enforcing uniqueness at EC level 3
     and excluding entries without accession identifiers (protein_refs).
  2) Finishes mapping by adding explicit H and assigning hydrogen atom maps using
     evodex.hydrogen_mapper.map_hydrogens_in_reaction.
  3) Extracts Bm, Cm, Em operators using evodex.operators.extract_operator_by_abstraction
     (with matched=True).
  4) Prints a per-enzyme summary of unique operator hashes and representative SMIRKS.
  5) Writes a per-enzyme CSV summary (counts per level) for manuscript text.

Outputs (in analysis/top5_diverse_enzymes_operators/output/):
  - top5_enzyme_mapped_reactions.csv
  - top5_enzyme_operators.csv
  - mechanistic_consistency_summary.csv
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

from evodex.hydrogen_mapper import map_hydrogens_in_reaction
from evodex.operators import extract_operator_by_abstraction
from evodex.utils import reaction_hash

EnzymeKey = Tuple[str, str, str, str]


def process_top5_enzyme_reactions(
    enzyme_reactions: List[Tuple[EnzymeKey, List[Dict[str, str]]]],
    output_file: str,
) -> None:
    """
    Add explicit H and assign hydrogen atom maps using hydrogen_mapper.

    Expects each input row to have a 'mapped' reaction string containing
    consistent heavy-atom maps on both sides of the reaction.
    """
    with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "organism",
            "protein_refs",
            "protein_db",
            "ec_num",
            "original_mapped",
            "remapped_smiles",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for enzyme_key, reactions in enzyme_reactions:
            organism, protein_refs, protein_db, ec_num = enzyme_key
            for row in reactions:
                try:
                    orig_smiles = row["mapped"]
                    remapped = map_hydrogens_in_reaction(orig_smiles)
                    writer.writerow(
                        {
                            "organism": organism,
                            "protein_refs": protein_refs,
                            "protein_db": protein_db,
                            "ec_num": ec_num,
                            "original_mapped": orig_smiles,
                            "remapped_smiles": remapped,
                        }
                    )
                except Exception as e:
                    print(f"Error processing reaction {row.get('rxn_idx', '')}: {e}")


def read_reaction_data(file_path: str) -> Dict[EnzymeKey, List[Dict[str, str]]]:
    enzyme_counts: Dict[EnzymeKey, List[Dict[str, str]]] = defaultdict(list)

    with open(file_path, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            enzyme_key: EnzymeKey = (
                row["organism"].strip(),
                row["protein_refs"].strip(),
                row["protein_db"].strip(),
                row["ec_num"].strip(),
            )
            enzyme_counts[enzyme_key].append(row)

    return enzyme_counts


def extract_bm_cm_em_operators(input_file: str, output_file: str) -> None:
    """
    Extract Bm, Cm, Em operators and hashes using extract_operator_by_abstraction.

    Mapping:
      - Bm = abstraction "B"
      - Cm = abstraction "C"
      - Em = abstraction "E"
    matched=True excludes unmapped atoms, matching prior "exclude unmapped" behavior.
    """
    levels = [("Bm", "B"), ("Cm", "C"), ("Em", "E")]

    with open(input_file, newline="", encoding="utf-8") as infile, open(
        output_file, "w", newline="", encoding="utf-8"
    ) as outfile:
        reader = csv.DictReader(infile)

        fieldnames = [
            "organism",
            "protein_refs",
            "protein_db",
            "ec_num",
            "original_mapped",
            "remapped_smiles",
            "Bm",
            "Bm_hash",
            "Cm",
            "Cm_hash",
            "Em",
            "Em_hash",
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            output_row = dict(row)
            rxn = row["remapped_smiles"]

            for label, level in levels:
                try:
                    op = extract_operator_by_abstraction(rxn, level, matched=True)
                    op_hash = reaction_hash(op)
                except Exception:
                    op = ""
                    op_hash = ""

                output_row[label] = op
                output_row[f"{label}_hash"] = op_hash

            writer.writerow(output_row)


def summarize_operator_hashes_by_enzyme(operator_file: str) -> None:
    enzyme_bins = defaultdict(lambda: {"Bm": Counter(), "Cm": Counter(), "Em": Counter()})
    row_lookup: Dict[Tuple[EnzymeKey, str, str], str] = {}

    with open(operator_file, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            enzyme_key: EnzymeKey = (
                row["organism"],
                row["protein_refs"],
                row["protein_db"],
                row["ec_num"],
            )

            for label in ["Bm", "Cm", "Em"]:
                h = row.get(f"{label}_hash", "")
                smirks = row.get(label, "")
                if h and smirks:
                    row_lookup[(enzyme_key, h, label)] = smirks
                    enzyme_bins[enzyme_key][label][h] += 1

    print("\nOperator Hash Summary by Enzyme\n" + "=" * 35)
    for idx, (enzyme_key, hashes) in enumerate(enzyme_bins.items(), 1):
        organism, protein_refs, protein_db, ec_num = enzyme_key
        print(f"\n[{idx}] {organism} — EC {ec_num}")
        if protein_refs or protein_db:
            print(f"     Protein: {protein_refs} @ {protein_db}")

        for label in ["Em", "Cm", "Bm"]:
            print(f"  {label} Operators ({len(hashes[label])} unique):")
            for h, count in hashes[label].most_common():
                smirks = row_lookup.get((enzyme_key, h, label), "")
                print(f"    • {h} (n={count})\n      {smirks}")


def write_mechanistic_consistency_summary(operator_file: str, output_file: str) -> None:
    """
    Write a per-enzyme summary CSV with counts needed to populate manuscript text.

    One row per enzyme, with:
      - n_reactions
      - number of unique operator hashes at Bm/Cm/Em
      - totals at each level (sanity check; typically equals n_reactions if extraction succeeds)
    """
    enzyme_bins = defaultdict(
        lambda: {
            "organism": None,
            "protein_refs": None,
            "protein_db": None,
            "ec_num": None,
            "n_reactions": 0,
            "Bm": Counter(),
            "Cm": Counter(),
            "Em": Counter(),
        }
    )

    with open(operator_file, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            key = (
                row["organism"],
                row["protein_refs"],
                row["protein_db"],
                row["ec_num"],
            )
            rec = enzyme_bins[key]

            rec["organism"] = row["organism"]
            rec["protein_refs"] = row["protein_refs"]
            rec["protein_db"] = row["protein_db"]
            rec["ec_num"] = row["ec_num"]
            rec["n_reactions"] += 1

            for level in ["Bm", "Cm", "Em"]:
                h = row.get(f"{level}_hash", "")
                if h:
                    rec[level][h] += 1

    with open(output_file, "w", newline="", encoding="utf-8") as out:
        fieldnames = [
            "organism",
            "ec_num",
            "protein_refs",
            "protein_db",
            "n_reactions",
            "Bm_unique",
            "Cm_unique",
            "Em_unique",
            "Bm_total",
            "Cm_total",
            "Em_total",
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()

        for rec in enzyme_bins.values():
            writer.writerow(
                {
                    "organism": rec["organism"],
                    "ec_num": rec["ec_num"],
                    "protein_refs": rec["protein_refs"],
                    "protein_db": rec["protein_db"],
                    "n_reactions": rec["n_reactions"],
                    "Bm_unique": len(rec["Bm"]),
                    "Cm_unique": len(rec["Cm"]),
                    "Em_unique": len(rec["Em"]),
                    "Bm_total": sum(rec["Bm"].values()),
                    "Cm_total": sum(rec["Cm"].values()),
                    "Em_total": sum(rec["Em"].values()),
                }
            )


def main() -> None:
    file_path = "data/raw/raw_reactions.csv"
    enzyme_counts = read_reaction_data(file_path)

    sorted_enzymes = sorted(enzyme_counts.items(), key=lambda item: len(item[1]), reverse=True)

    print("Top 5 most diverse enzymes by unique EC level 3:\n")
    seen_ec_level3 = set()
    diverse_enzymes: List[Tuple[EnzymeKey, List[Dict[str, str]]]] = []

    for enzyme_key, reactions in sorted_enzymes:
        ec_full = enzyme_key[3]
        ec_level3 = ".".join(ec_full.split(".")[:3]) if ec_full else ""
        protein_refs = enzyme_key[1]

        if ec_level3 and (ec_level3 not in seen_ec_level3) and protein_refs.strip("[]").strip():
            seen_ec_level3.add(ec_level3)
            diverse_enzymes.append((enzyme_key, reactions))

        if len(diverse_enzymes) == 5:
            break

    for i, (enzyme_key, reactions) in enumerate(diverse_enzymes, 1):
        organism, protein_refs, protein_db, ec_num = enzyme_key
        print(
            f"{i}. Organism: {organism}, Protein Refs: {protein_refs}, "
            f"DB: {protein_db}, EC: {ec_num} — {len(reactions)} reactions"
        )

    mapped_out = "analysis/top5_diverse_enzymes_operators/output/top5_enzyme_mapped_reactions.csv"
    process_top5_enzyme_reactions(diverse_enzymes, mapped_out)
    print(f"\nRemapped reaction SMILES written to {mapped_out}")

    ops_out = "analysis/top5_diverse_enzymes_operators/output/top5_enzyme_operators.csv"
    extract_bm_cm_em_operators(mapped_out, ops_out)
    print(f"Operator extraction complete. Results written to {ops_out}")

    summary_out = "analysis/top5_diverse_enzymes_operators/output/mechanistic_consistency_summary.csv"
    write_mechanistic_consistency_summary(ops_out, summary_out)
    print(f"Mechanistic consistency summary written to {summary_out}")

    summarize_operator_hashes_by_enzyme(ops_out)


if __name__ == "__main__":
    main()
