#!/usr/bin/env python3
"""
TwinMap: deterministic atom mapping transfer between substrate and product using an operator.

Inputs
  - Substrate: SMILES
  - Product: SMILES
  - Operator: SMIRKS with optional atom maps

Output
  A JSON object with:
    - mapped_atoms: (([(idx, map)], [(idx, map)])) for substrate and product
    - unmapped_matched_atoms: (([idx], [idx])) for substrate and product
    - unmatched_atoms: (([idx], [idx])) for substrate and product
    - pairs: [(substrate_idx, product_idx, map)]
    - mapped_substrate_smiles: SMILES with atom map numbers
    - mapped_product_smiles: SMILES with atom map numbers

Usage
  python twinmap.py --substrate "CCO" --product "CC=O" --operator "[C:1][O:2]>>[C:1]=[O:2]"
"""

import json
import argparse
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import rdChemReactions
try:
    from rdkit.Chem import rdMolDraw2D  # optional (Cairo)
    HAVE_RDMOLDRAW2D = True
except Exception:
    rdMolDraw2D = None
    HAVE_RDMOLDRAW2D = False
from rdkit.Chem import rdDepictor, Draw
from collections import deque, defaultdict


def _reduce_mol(mol: Chem.Mol, remove_indices: Tuple[int, ...]) -> str:
    rw = Chem.RWMol()
    idx_map = {}
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        if i not in remove_indices:
            idx_map[i] = rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in idx_map and b in idx_map:
            rw.AddBond(idx_map[a], idx_map[b], bond.GetBondType())
    return Chem.MolToSmiles(rw.GetMol(), canonical=True)


def _mol_from_smiles(smiles: str) -> Chem.Mol:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"Bad SMILES: {smiles}")
    Chem.SanitizeMol(m)
    return m


def _smirks_to_rxn(smirks: str) -> rdChemReactions.ChemicalReaction:
    rxn = rdChemReactions.ReactionFromSmarts(smirks, useSmiles=False)
    if rxn is None:
        raise ValueError(f"Bad SMIRKS: {smirks}")
    if rxn.GetNumReactantTemplates() != 1 or rxn.GetNumProductTemplates() != 1:
        raise ValueError("TwinMap expects a unary operator with 1 reactant and 1 product template")
    return rxn


def _draw_with_maps(smiles: str, out_png: str, width: int = 600, height: int = 450) -> None:
    """Draw a molecule PNG labeling atoms by their atom-map numbers only.
    Uses rdMolDraw2D if available, otherwise falls back to Chem.Draw (Pillow).
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"Cannot draw; bad SMILES: {smiles}")
    m = Chem.RemoveHs(m)
    rdDepictor.Compute2DCoords(m)

    # Map-number-only labels
    atom_labels = {}
    for a in m.GetAtoms():
        amap = a.GetAtomMapNum()
        atom_labels[a.GetIdx()] = (str(amap) if amap > 0 else "")

    if HAVE_RDMOLDRAW2D:
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        opts = drawer.drawOptions()
        opts.atomLabels = atom_labels
        drawer.DrawMolecule(m)
        drawer.FinishDrawing()
        png = drawer.GetDrawingText()
        with open(out_png, "wb") as f:
            f.write(png)
        return

    # Fallback: use Chem.Draw (Pillow). Encode labels via the 'atomLabel' atom property.
    m2 = Chem.Mol(m)
    for idx, label in atom_labels.items():
        atom = m2.GetAtomWithIdx(idx)
        if label:
            atom.SetProp('atomLabel', label)
        else:
            try:
                atom.ClearProp('atomLabel')
            except Exception:
                pass
    img = Draw.MolToImage(m2, size=(width, height))
    img.save(out_png)


def _assign_seq_maps_on_unmapped_pairs(
    rxn_reactant: Chem.Mol,
    rxn_product: Chem.Mol,
    op_reactant: Chem.Mol,
    op_product: Chem.Mol,
    final_r_match: Tuple[int, ...],
    final_p_match: Tuple[int, ...],
) -> None:
    # Copy operator maps to matched atoms that already have maps in the operator
    for i, r_idx in enumerate(final_r_match):
        mnum = op_reactant.GetAtomWithIdx(i).GetAtomMapNum()
        if mnum > 0:
            rxn_reactant.GetAtomWithIdx(r_idx).SetAtomMapNum(mnum)
    for i, p_idx in enumerate(final_p_match):
        mnum = op_product.GetAtomWithIdx(i).GetAtomMapNum()
        if mnum > 0:
            rxn_product.GetAtomWithIdx(p_idx).SetAtomMapNum(mnum)

    # Determine next map number
    existing = []
    for a in op_reactant.GetAtoms():
        if a.GetAtomMapNum() > 0:
            existing.append(a.GetAtomMapNum())
    for a in op_product.GetAtoms():
        if a.GetAtomMapNum() > 0:
            existing.append(a.GetAtomMapNum())
    next_map = (max(existing) if existing else 0) + 1

    # Assign new maps to positions where both sides are matched but unmapped in the operator
    for i, (r_idx, p_idx) in enumerate(zip(final_r_match, final_p_match)):
        lhs_map = op_reactant.GetAtomWithIdx(i).GetAtomMapNum()
        rhs_map = op_product.GetAtomWithIdx(i).GetAtomMapNum()
        if lhs_map == 0 and rhs_map == 0:
            if rxn_reactant.GetAtomWithIdx(r_idx).GetAtomMapNum() == 0 and \
               rxn_product.GetAtomWithIdx(p_idx).GetAtomMapNum() == 0:
                rxn_reactant.GetAtomWithIdx(r_idx).SetAtomMapNum(next_map)
                rxn_product.GetAtomWithIdx(p_idx).SetAtomMapNum(next_map)
                next_map += 1


# ------- Extend maps to the unchanged region via synchronized BFS -------

essential_bond_types = {
    Chem.BondType.SINGLE,
    Chem.BondType.DOUBLE,
    Chem.BondType.TRIPLE,
    Chem.BondType.AROMATIC,
}


def _atom_sig(m: Chem.Mol, idx: int) -> tuple:
    """Local, deterministic signature used to sort candidates and break symmetries."""
    a = m.GetAtomWithIdx(idx)
    nbr_atomic = sorted(n.GetAtomicNum() for n in a.GetNeighbors())
    return (
        a.GetAtomicNum(),
        a.GetFormalCharge(),
        int(a.GetIsAromatic()),
        a.GetTotalDegree(),
        tuple(nbr_atomic),
        idx,  # tie-breaker for determinism
    )


def _extend_maps_by_bfs_propagation(rxn_reactant: Chem.Mol, rxn_product: Chem.Mol) -> None:
    """Extend atom-map numbers from the reaction center to the rest of the molecule.

    Seeds are atoms that already share the same map number. For each seed pair (r_i, p_i),
    pair *unmapped* neighbors connected by the same bond order and with the same element
    and aromaticity. Pairing is made deterministic by sorting by a local signature.
    Assign fresh map numbers to the paired neighbors and continue BFS until convergence.
    """
    # Build seed pairs from atoms that already share a map on both sides
    by_map_r = {}
    by_map_p = {}
    for a in rxn_reactant.GetAtoms():
        mnum = a.GetAtomMapNum()
        if mnum > 0:
            by_map_r[mnum] = a.GetIdx()
    for a in rxn_product.GetAtoms():
        mnum = a.GetAtomMapNum()
        if mnum > 0:
            by_map_p[mnum] = a.GetIdx()

    seeds = [(r_idx, by_map_p[m]) for m, r_idx in by_map_r.items() if m in by_map_p]
    if not seeds:
        return

    # Next map number after any existing assignments on either side
    existing = [a.GetAtomMapNum() for a in rxn_reactant.GetAtoms() if a.GetAtomMapNum() > 0]
    existing += [a.GetAtomMapNum() for a in rxn_product.GetAtoms() if a.GetAtomMapNum() > 0]
    next_map = (max(existing) if existing else 0) + 1

    q = deque(seeds)
    while q:
        r_i, p_i = q.popleft()

        # Collect *unmapped* neighbors of r_i and p_i (grouped by bond type)
        r_by_bt: defaultdict = defaultdict(list)
        for b in rxn_reactant.GetAtomWithIdx(r_i).GetBonds():
            bt = b.GetBondType()
            if bt not in essential_bond_types:
                continue
            j = b.GetOtherAtomIdx(r_i)
            if rxn_reactant.GetAtomWithIdx(j).GetAtomMapNum() == 0:
                r_by_bt[bt].append(j)

        p_by_bt: defaultdict = defaultdict(list)
        for b in rxn_product.GetAtomWithIdx(p_i).GetBonds():
            bt = b.GetBondType()
            if bt not in essential_bond_types:
                continue
            j = b.GetOtherAtomIdx(p_i)
            if rxn_product.GetAtomWithIdx(j).GetAtomMapNum() == 0:
                p_by_bt[bt].append(j)

        # For each shared bond type, pair neighbors deterministically and assign maps
        for bt in list(set(r_by_bt.keys()) & set(p_by_bt.keys())):
            r_list = sorted(r_by_bt[bt], key=lambda j: _atom_sig(rxn_reactant, j))
            p_list = sorted(p_by_bt[bt], key=lambda j: _atom_sig(rxn_product, j))
            for r_j, p_j in zip(r_list, p_list):
                ar = rxn_reactant.GetAtomWithIdx(r_j)
                ap = rxn_product.GetAtomWithIdx(p_j)
                if ar.GetAtomicNum() != ap.GetAtomicNum():
                    continue
                if ar.GetIsAromatic() != ap.GetIsAromatic():
                    continue
                if ar.GetAtomMapNum() == 0 and ap.GetAtomMapNum() == 0:
                    ar.SetAtomMapNum(next_map)
                    ap.SetAtomMapNum(next_map)
                    q.append((r_j, p_j))
                    next_map += 1


def _collect_bins(
    rxn_reactant: Chem.Mol,
    rxn_product: Chem.Mol,
    op_reactant: Chem.Mol,
    op_product: Chem.Mol,
    final_r_match: Tuple[int, ...],
    final_p_match: Tuple[int, ...],
):
    mapped_r, mapped_p = [], []
    unmapped_r, unmapped_p = [], []

    for i, r_idx in enumerate(final_r_match):
        mnum = op_reactant.GetAtomWithIdx(i).GetAtomMapNum()
        if mnum > 0:
            mapped_r.append((r_idx, mnum))
        else:
            mnum_rxn = rxn_reactant.GetAtomWithIdx(r_idx).GetAtomMapNum()
            if mnum_rxn > 0:
                mapped_r.append((r_idx, mnum_rxn))
            else:
                unmapped_r.append(r_idx)

    for i, p_idx in enumerate(final_p_match):
        mnum = op_product.GetAtomWithIdx(i).GetAtomMapNum()
        if mnum > 0:
            mapped_p.append((p_idx, mnum))
        else:
            mnum_rxn = rxn_product.GetAtomWithIdx(p_idx).GetAtomMapNum()
            if mnum_rxn > 0:
                mapped_p.append((p_idx, mnum_rxn))
            else:
                unmapped_p.append(p_idx)

    all_r = set(range(rxn_reactant.GetNumAtoms()))
    all_p = set(range(rxn_product.GetNumAtoms()))
    matched_r = set([i for i, _ in mapped_r]) | set(unmapped_r)
    matched_p = set([i for i, _ in mapped_p]) | set(unmapped_p)
    unmatched = (sorted(all_r - matched_r), sorted(all_p - matched_p))

    by_map_r = {m: i for i, m in mapped_r}
    pairs = []
    used_p = set()
    for i, m in mapped_r:
        if m in by_map_r:
            j = None
            for ip, mp in mapped_p:
                if mp == m:
                    j = ip
                    break
            if j is not None:
                pairs.append((i, j, m))
                used_p.add(j)

    # Note: newly BFS-mapped atoms are not in final_*_match, but they are present in the
    # SMILES we emit below; update of `pairs` beyond the match is left minimal by design.

    return (mapped_r, mapped_p), (unmapped_r, unmapped_p), unmatched, pairs


def twinmap(
    substrate_smiles: str,
    product_smiles: str,
    operator_smirks: str,
) -> Dict:
    # Build reaction with one reactant and one product using the provided SMILES
    reactant = _mol_from_smiles(substrate_smiles)
    product = _mol_from_smiles(product_smiles)

    rxn = rdChemReactions.ChemicalReaction()
    rxn.AddReactantTemplate(Chem.Mol(reactant))
    rxn.AddProductTemplate(Chem.Mol(product))

    op = _smirks_to_rxn(operator_smirks)
    op_r = op.GetReactantTemplate(0)
    op_p = op.GetProductTemplate(0)

    # Find matches
    r_matches = reactant.GetSubstructMatches(op_r, uniquify=False)
    p_matches = product.GetSubstructMatches(op_p, uniquify=False)

    # Reduced forms for pairing
    reduced_r = [_reduce_mol(reactant, m) for m in r_matches]
    reduced_p = [_reduce_mol(product, m) for m in p_matches]

    final_r = None
    final_p = None
    for i, rs in enumerate(reduced_r):
        for j, ps in enumerate(reduced_p):
            if rs == ps:
                final_r = r_matches[i]
                final_p = p_matches[j]
                break
        if final_r is not None:
            break

    if final_r is None or final_p is None:
        return {
            "mapped_atoms": ([], []),
            "unmapped_matched_atoms": ([], []),
            "unmatched_atoms": (
                list(range(reactant.GetNumAtoms())),
                list(range(product.GetNumAtoms()))
            ),
            "pairs": [],
            "mapped_substrate_smiles": Chem.MolToSmiles(Chem.RemoveHs(reactant)),
            "mapped_product_smiles": Chem.MolToSmiles(Chem.RemoveHs(product)),
        }

    # Work on copies so we can set map numbers
    reactant_c = Chem.Mol(reactant)
    product_c = Chem.Mol(product)

    _assign_seq_maps_on_unmapped_pairs(
        reactant_c, product_c, op_r, op_p, final_r, final_p
    )

    # NEW: extend maps to the unchanged region
    _extend_maps_by_bfs_propagation(reactant_c, product_c)

    mapped, unmapped, unmatched, pairs = _collect_bins(
        reactant_c, product_c, op_r, op_p, final_r, final_p
    )

    # Final SMILES with map numbers
    r_noh = Chem.RemoveHs(reactant_c)
    p_noh = Chem.RemoveHs(product_c)
    mapped_substrate_smiles = Chem.MolToSmiles(r_noh)
    mapped_product_smiles = Chem.MolToSmiles(p_noh)

    return {
        "mapped_atoms": mapped,
        "unmapped_matched_atoms": unmapped,
        "unmatched_atoms": unmatched,
        "pairs": pairs,
        "mapped_substrate_smiles": mapped_substrate_smiles,
        "mapped_product_smiles": mapped_product_smiles,
    }


def main():
    ap = argparse.ArgumentParser(description="TwinMap: map substrate and product using an operator")
    ap.add_argument("--substrate", required=True, help="Substrate SMILES")
    ap.add_argument("--product", required=True, help="Product SMILES")
    ap.add_argument("--operator", required=True, help="Operator SMIRKS")
    ap.add_argument("--pretty", action="store_true", help="Pretty print JSON")
    ap.add_argument("--png-prefix", help="If set, write PNGs of mapped substrate and product using this prefix")
    args = ap.parse_args()

    res = twinmap(args.substrate, args.product, args.operator)
    if args.pretty:
        print(json.dumps(res, indent=2))
    else:
        print(json.dumps(res, separators=(",", ":")))

    if args.png_prefix:
        try:
            sub_png = f"{args.png_prefix}_substrate.png"
            prod_png = f"{args.png_prefix}_product.png"
            _draw_with_maps(res["mapped_substrate_smiles"], sub_png)
            _draw_with_maps(res["mapped_product_smiles"], prod_png)
            print(f"\nWrote: {sub_png}\nWrote: {prod_png}")
        except Exception as e:
            print(f"\n[warn] Failed to render PNGs: {e}")


if __name__ == "__main__":
    main()

#  python twinmap.py \
#   --substrate "OCC(C)CCCCO" \
#   --product   "OCC(C)CCCC=O" \
#   --operator  "[C:1][O:2]>>[C:1]=[O:2]" \
#   --pretty \
#   --png-prefix example