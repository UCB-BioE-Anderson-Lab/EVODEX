#!/usr/bin/env python3
"""
Convert EVODEX CSV tables into JSON for the website.

- Inputs (hard-coded relative paths):
    EVODEX/evodex/data/*.csv

- Outputs:
    EVODEX/website/data/<CSV_STEM>.json

Behavior:
- For any CSV that has a 'smirks' column (A, Am, ..., E, Em, P, R, E subsets, etc.),
  each JSON entry is keyed by the CSV 'id' field and contains:
    {
      "<EVODEX-ID>": {
        "smirks": "...original SMIRKS...",
        "substrates": [... isotope-encoded SMILES ...],
        "products":  [... isotope-encoded SMILES ...],
        "sources":   [... parsed source IDs ...],
        "meta":      {<any additional columns from the CSV, if present>}
      }
    }

- For CSVs without a 'smirks' column (e.g. EVODEX-M, EVODEX-F, selected_reactions),
  each JSON entry is keyed by 'id' and the value is simply the remaining CSV fields
  from that row (no special processing).
"""

import csv
import json
import os
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem import AllChem


def mol_to_isotope_smiles(mol: Chem.Mol) -> str:
    """
    Copy a template mol and convert atom-map numbers into isotopes so that
    SmilesDrawer will render map indices as isotope labels.

    [C:23] -> [23C]
    """
    mol = Chem.Mol(mol)  # shallow copy

    for atom in mol.GetAtoms():
        # Move atom-map numbers into isotopes (so SmilesDrawer shows them as labels)
        amap = atom.GetAtomMapNum()
        if amap:
            atom.SetIsotope(amap)
            atom.SetAtomMapNum(0)

        # For non-hydrogen atoms, suppress implicit/aggregated hydrogens.
        # We only want hydrogens that were explicitly present in the SMIRKS
        # (as separate [#1:map] atoms) to appear. Otherwise RDKit will
        # generate things like [23CH3], which show up as "H3" on the carbon.
        if atom.GetAtomicNum() != 1:
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)

    # Keep stereochemistry if present
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def split_smirks_to_smiles(smirks: str):
    """
    Given a SMIRKS string, return (substrate_smiles_list, product_smiles_list).

    Uses RDKit's reaction SMARTS parser, then converts each reactant and
    product template to a SMILES where the atom-map numbers are encoded
    as isotopes.
    """
    try:
        rxn = AllChem.ReactionFromSmarts(smirks, useSmiles=False)
    except Exception as e:
        raise ValueError(f"Failed to parse SMIRKS: {smirks}") from e

    substrates = []
    products = []

    for i in range(rxn.GetNumReactantTemplates()):
        tmpl = rxn.GetReactantTemplate(i)
        substrates.append(mol_to_isotope_smiles(tmpl))

    for i in range(rxn.GetNumProductTemplates()):
        tmpl = rxn.GetProductTemplate(i)
        products.append(mol_to_isotope_smiles(tmpl))

    return substrates, products


def parse_sources_field(raw) -> list[str]:
    """
    Turn the CSV 'sources' field into a list of IDs.

    Handles single IDs or multiple IDs separated by ',' or ';'.
    """
    if raw is None:
        return []
    s = str(raw).strip()
    if not s:
        return []
    parts = re.split(r"[;,]", s)
    return [p.strip() for p in parts if p.strip()]


def build_paths():
    """
    Compute input/output base directories relative to this script.

    Script is at:
        EVODEX/pipeline/phase7_website.py

    We want:
        input dir:  EVODEX/evodex/data/
        output dir: EVODEX/website/data/
    """
    script_path = Path(__file__).resolve()
    repo_root = script_path.parent.parent  # .. from pipeline/

    evodex_data_dir = repo_root / "evodex" / "data"
    website_data_dir = repo_root / "website" / "data"

    return evodex_data_dir, website_data_dir


def csv_to_json():
    """
    Convert all EVODEX CSVs under evodex/data into JSON files under website/data.
    """
    evodex_data_dir, website_data_dir = build_paths()

    if not evodex_data_dir.is_dir():
        raise FileNotFoundError(f"EVODEX data directory not found: {evodex_data_dir}")

    website_data_dir.mkdir(parents=True, exist_ok=True)

    for csv_path in sorted(evodex_data_dir.glob("*.csv")):
        output: dict[str, dict] = {}

        with csv_path.open(newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames or []
            has_smirks = "smirks" in fieldnames

            for row in reader:
                eid = str(row.get("id", "")).strip()
                if not eid:
                    continue

                if has_smirks:
                    smirks = (row.get("smirks") or "").strip()
                    if not smirks:
                        continue

                    try:
                        substrates, products = split_smirks_to_smiles(smirks)
                    except Exception as e:
                        print(f"[WARN] {csv_path.name}: skipping {eid}: {e}")
                        continue

                    sources_raw = row.get("sources", "")
                    entry = {
                        "smirks": smirks,
                        "substrates": substrates,
                        "products": products,
                        "sources": parse_sources_field(sources_raw),
                    }

                    # Preserve any additional columns as metadata
                    meta = {
                        k: v
                        for k, v in row.items()
                        if k not in ("id", "smirks", "sources") and v not in (None, "")
                    }
                    if meta:
                        entry["meta"] = meta

                    output[eid] = entry
                else:
                    # Generic CSV: just map id -> remaining fields
                    value = {k: v for k, v in row.items() if k != "id"}
                    output[eid] = value

        json_path = website_data_dir / f"{csv_path.stem}.json"
        with json_path.open("w") as jf:
            json.dump(output, jf, indent=2)

        print(f"Wrote {len(output)} entries from {csv_path.name} to {json_path}")


if __name__ == "__main__":
    csv_to_json()