#!/usr/bin/env python3

from pathlib import Path
import json

from rdkit import Chem
from rdkit.Chem import AllChem

from evodex.utils import reaction_hash


BASE_DIR = Path(__file__).resolve().parent
IN_DIR = BASE_DIR / "out/2_map"
OUT_DIR = BASE_DIR / "out/3_hmap"
LOG_DIR = BASE_DIR / "out/_logs"

# If non-empty, restrict processing to these table_ids
ONLY_TABLES = []  # e.g. ["trifluoroacetanilide_methanolysis"]


def stream_jsonl(path: Path):
    """Yield one JSON object per non-empty line."""
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)


def mol_from_smiles_checked(smiles: str) -> Chem.Mol:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"Bad SMILES: {smiles}")
    Chem.SanitizeMol(m)
    return m


def add_hydrogens_and_convert_to_at(smiles: str) -> str:
    """Add explicit H and then convert all H (atomic number 1) to At (85).

    Keeps existing atom-map numbers on heavy atoms.
    """
    m = mol_from_smiles_checked(smiles)
    m_h = Chem.AddHs(m)
    rw = Chem.RWMol(m_h)
    for a in rw.GetAtoms():
        if a.GetAtomicNum() == 1:  # hydrogen
            a.SetAtomicNum(85)  # astatine
    return Chem.MolToSmiles(rw.GetMol(), kekuleSmiles=False)


def convert_at_to_h(smiles: str) -> str:
    """Convert At (85) back to H (1), preserving atom-map numbers."""
    m = mol_from_smiles_checked(smiles)
    rw = Chem.RWMol(m)
    for a in rw.GetAtoms():
        if a.GetAtomicNum() == 85:
            a.SetAtomicNum(1)
    return Chem.MolToSmiles(rw.GetMol(), kekuleSmiles=False)


# Relaxed mapping: assign At maps based on heavy-atom maps, without strict checking
def map_atoms_relaxed(smirks: str) -> str:
    """
    Assign atom maps to astatines (representing hydrogens) based on existing
    heavy-atom maps, without enforcing that *all* heavy atoms are mapped or
    that reactant/product map sets are identical.

    This is a relaxed variant of the EVODEX map_atoms logic that:
      - assumes heavy-atom maps already exist on the reaction core
      - uses those maps to deterministically assign H/At maps
      - leaves heavy-atom maps untouched
    """
    rxn = AllChem.ReactionFromSmarts(smirks, useSmiles=True)
    if rxn is None:
        raise ValueError(f"Invalid reaction SMIRKS: {smirks}")

    # Work on local copies of reactants and products
    reactants = [Chem.Mol(r) for r in rxn.GetReactants()]
    products = [Chem.Mol(p) for p in rxn.GetProducts()]

    # Collect existing heavy-atom map numbers
    all_maps = set()
    for mol in reactants + products:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 85:  # At
                continue
            mnum = atom.GetAtomMapNum()
            if mnum > 0:
                all_maps.add(mnum)

    if not all_maps:
        raise ValueError(f"No existing heavy-atom atom maps found in {smirks}")

    next_map = max(all_maps) + 1

    # For each heavy-atom map, we will keep a queue of H/At map numbers
    # assigned on the reactant side, to reuse on the product side.
    astatine_map_dict = {}  # heavy_map -> [hmap1, hmap2, ...]

    def _add_astatine_maps(mol: Chem.Mol, is_reactant: bool) -> None:
        nonlocal next_map
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 85:
                continue
            # Attach At to the first non-At neighbor with a valid map, if any
            neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 85]
            if not neighbors:
                continue
            neighbor_map = neighbors[0].GetAtomMapNum()
            if neighbor_map <= 0:
                # Neighbor has no map; skip assigning a map to this At
                continue

            if is_reactant:
                astatine_map_dict.setdefault(neighbor_map, []).append(next_map)
                atom.SetAtomMapNum(next_map)
                next_map += 1
            else:
                lst = astatine_map_dict.get(neighbor_map)
                if lst:
                    atom.SetAtomMapNum(lst.pop(0))
                else:
                    # No corresponding H on the reactant side; leave unmapped
                    atom.SetAtomMapNum(0)

    for r in reactants:
        _add_astatine_maps(r, is_reactant=True)
    for p in products:
        _add_astatine_maps(p, is_reactant=False)

    # Convert back to SMIRKS string
    reactant_smarts = [Chem.MolToSmarts(m, isomericSmiles=True) for m in reactants]
    product_smarts = [Chem.MolToSmarts(m, isomericSmiles=True) for m in products]
    return ">>".join([".".join(reactant_smarts), ".".join(product_smarts)])


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    errors_path = LOG_DIR / "3_hmap.errors.jsonl"
    with errors_path.open("w", encoding="utf-8") as err_f:
        input_files = sorted(IN_DIR.glob("*.atommap.jsonl"))

        for in_path in input_files:
            table_id = in_path.name.replace(".atommap.jsonl", "")

            if ONLY_TABLES and table_id not in ONLY_TABLES:
                continue

            out_path = OUT_DIR / f"{table_id}.hmap.jsonl"
            with out_path.open("w", encoding="utf-8") as out_f:
                for row in stream_jsonl(in_path):
                    sub_mapped = row.get("mapped_substrate_smiles")
                    prod_mapped = row.get("mapped_product_smiles")

                    if not sub_mapped or not prod_mapped:
                        err_rec = {
                            "table_id": table_id,
                            "row": row.get("row"),
                            "label": row.get("label"),
                            "stage": "input",
                            "error": "missing_mapped_substrate_or_product",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    try:
                        # 1) Add explicit H and convert H -> At for both sides
                        sub_at = add_hydrogens_and_convert_to_at(sub_mapped)
                        prod_at = add_hydrogens_and_convert_to_at(prod_mapped)

                        # 2) Assign maps to At based on heavy-atom maps (relaxed mapping)
                        smirks_at = f"{sub_at}>>{prod_at}"
                        mapped_smirks_at = map_atoms_relaxed(smirks_at)

                        # 3) Convert At back to H in the mapped SMIRKS
                        try:
                            left_at, right_at = mapped_smirks_at.split(">>", 1)
                        except ValueError:
                            raise ValueError(f"Could not split mapped SMIRKS: {mapped_smirks_at}")

                        sub_h = convert_at_to_h(left_at)
                        prod_h = convert_at_to_h(right_at)

                        mapped_smirks = f"{sub_h}>>{prod_h}"
                        rxn_hash = reaction_hash(mapped_smirks)

                    except Exception as e:
                        err_rec = {
                            "table_id": table_id,
                            "row": row.get("row"),
                            "label": row.get("label"),
                            "stage": "hydrogen_mapping",
                            "error": f"{e}",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    # Build a minimal Stage-3 record: keep identifiers, per-row data,
                    # and hydrogen-mapped EVODEX-P fields, but drop bulky intermediates.
                    rec = {
                        "table_id": row.get("table_id", table_id),
                        "row": row.get("row"),
                        "label": row.get("label"),
                        "data": row.get("data"),
                        # Stage-2 SMILES context
                        "substrate_smiles": row.get("substrate_smiles"),
                        "product_smiles": row.get("product_smiles"),
                        "mapped_substrate_smiles": sub_mapped,
                        "mapped_product_smiles": prod_mapped,
                        # New Stage-3 EVODEX-P hydrogen-mapped fields
                        "evodex_p_smirks": mapped_smirks,
                        "evodex_p_hash": rxn_hash,
                        "evodex_p_substrate_smiles": sub_h,
                        "evodex_p_product_smiles": prod_h,
                    }

                    out_f.write(json.dumps(rec, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
