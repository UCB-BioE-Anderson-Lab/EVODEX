from typing import Dict, Tuple, List, Optional, Any
from collections import Counter
from pathlib import Path
import csv

from rdkit import Chem
from rdkit.Chem import rdChemReactions, AllChem

# Allow very large CSV fields (needed for long SMIRKS strings)
csv.field_size_limit(10**7)


def prepare_operator(operator_smirks: str):
    """
    Build an RDKit ChemicalReaction object from an operator SMIRKS or SMARTS.

    Parameters
    ----------
    operator_smirks : str
        SMIRKS or SMARTS string defining the reaction operator.

    Returns
    -------
    rdkit.Chem.rdChemReactions.ChemicalReaction
        Parsed operator as an RDKit ChemicalReaction.
    """
    operator = rdChemReactions.ReactionFromSmarts(operator_smirks)
    return operator


def prepare_reaction(reaction_smiles: str):
    """
    Prepare an RDKit ChemicalReaction object from a reaction SMILES string.

    This function only supports single-molecule reactions (one substrate to one product).
    If the reaction SMILES contains multiple molecules (indicated by '.' separators),
    a ValueError will be raised.

    Parameters
    ----------
    reaction_smiles : str
        A reaction SMILES string of the form 'reactant>>product' where both
        reactant and product must be single molecules.

    Returns
    -------
    rdkit.Chem.rdChemReactions.ChemicalReaction
        An RDKit ChemicalReaction object with hydrogens added to the
        reactant and product molecules.

    Raises
    ------
    ValueError
        If the reaction SMILES contains '.' indicating multiple molecules
        on either side of the reaction.
    """
    if '.' in reaction_smiles:
        raise ValueError(
            "Multi-molecule reactions are not supported. "
            "Reaction SMILES must be single substrate to single product. "
            f"Got: {reaction_smiles}"
        )
    
    reactant_smiles, product_smiles = reaction_smiles.split(">>")
    
    reactant_mol = Chem.AddHs(Chem.MolFromSmiles(reactant_smiles))
    product_mol = Chem.AddHs(Chem.MolFromSmiles(product_smiles))

    rxn = AllChem.ChemicalReaction()
    rxn.AddReactantTemplate(reactant_mol)
    rxn.AddProductTemplate(product_mol)

    return rxn


def mechanistic_label_reaction(rxn, operator) -> Dict:
    """
    Perform a relaxed, mechanistically oriented alignment between a concrete reaction
    and a reaction operator.

    This function treats the reaction as a single-substrate and single-product
    transformation and uses the first reactant and first product templates in each
    ChemicalReaction. It identifies atom-level correspondences between the operator
    and the reaction and classifies reaction atoms by their parity in the transformation.

    In contrast to a strict complete-operator reaction_labeler that requires the
    non-transformed remainder of the substrate and product to be identical, this
    mechanistic labeler allows extra material on either side. Thus it can be used
    on both matched and complete forms of partial reaction operators, such both
    E and Em. 
    
    The key steps are:

      1. Substructure matching:
         - Find all substructure matches of the operator reactant in the reaction
           substrate.
         - Find all substructure matches of the operator product in the reaction
           product.
         - If either side has no matches, the function returns an empty labeling.

      2. Fragment-based pairing:
         - For each substrate match, remove the matched atoms and decompose the
           remaining graph into disconnected fragments. Represent each fragment
           as a canonical SMILES.
         - Do the same for each product match.
         - For every combination of one substrate match and one product match,
           compare the multisets of fragment SMILES. A pairing is accepted if
           the fragments on one side form a multiset subset of the fragments on
           the other side (that is, p subset r or r subset p). This subset
           criterion allows the operator to be applied in contexts where
           additional groups are present on only one side of the reaction.

      3. Atom classification:
         - For the accepted pairing (if any), transfer atom-map numbers from the
           operator onto the corresponding atoms in the reaction:
             * mapped atoms: reaction atoms whose operator counterparts carry a
               positive atom map number;
             * unmapped matched atoms: reaction atoms that are structurally
               matched by the operator but lack an atom map number in the
               operator;
             * unmatched atoms: all remaining atoms in the substrate and product.

    Parameters
    ----------
    rxn : rdChemReactions.ChemicalReaction
        Reaction that may contain one or more reactant and product templates.
        The atom indices in the returned labels refer to
        rxn.GetReactants()[0] and rxn.GetProducts()[0].
    operator : rdChemReactions.ChemicalReaction
        Reaction operator that may contain one or more reactant and product
        templates. This function uses operator.GetReactants()[0] and
        operator.GetProducts()[0] when computing labels.

    Returns
    -------
    Dict
        A dictionary with the following keys:

        - "rxn": the original rxn object, which serves as the reference for
          atom indices.
        - "mapped_atoms": tuple of two lists, (reactant_mapped, product_mapped),
          where each list contains (atom_idx, map_num) pairs.
        - "unmapped_matched_atoms": tuple of two lists,
          (reactant_unmapped, product_unmapped), giving atom indices that are
          matched by the operator but have no map number.
        - "unmatched_atoms": tuple of two lists,
          (reactant_unmatched, product_unmatched), giving atom indices not
          involved in the operator match.

        If no acceptable substrate to product pairing is found, all three sets of
        atom lists are returned empty.
    """

    rxn_reactant = rxn.GetReactants()[0]
    rxn_product = rxn.GetProducts()[0]

    operator_reactant = operator.GetReactants()[0]
    operator_product = operator.GetProducts()[0]

    # All substructure matches of operator in substrate and product
    reactant_matches = rxn_reactant.GetSubstructMatches(
        operator_reactant, uniquify=False
    )
    product_matches = rxn_product.GetSubstructMatches(
        operator_product, uniquify=False
    )

    # Default empty labels that can be returned directly in early exit cases
    mapped_atoms: Tuple[List[Tuple[int, int]], List[Tuple[int, int]]] = ([], [])
    unmapped_matched_atoms: Tuple[List[int], List[int]] = ([], [])
    unmatched_atoms: Tuple[List[int], List[int]] = ([], [])

    # If either side has no matches, return an empty labeling
    if not reactant_matches or not product_matches:
        return {
            "rxn": rxn,
            "mapped_atoms": mapped_atoms,
            "unmapped_matched_atoms": unmapped_matched_atoms,
            "unmatched_atoms": unmatched_atoms,
        }

    # Precompute fragment lists for each match
    reduced_reactant_frags: List[List[str]] = []
    for reactant_match in reactant_matches:
        frags = _reduced_fragment_smiles(rxn_reactant, reactant_match)
        reduced_reactant_frags.append(frags)

    reduced_product_frags: List[List[str]] = []
    for product_match in product_matches:
        frags = _reduced_fragment_smiles(rxn_product, product_match)
        reduced_product_frags.append(frags)

    # Find best pairing using subset-of-fragments criterion
    final_reactant_match = None
    final_product_match = None

    for i, r_frags in enumerate(reduced_reactant_frags):
        for j, p_frags in enumerate(reduced_product_frags):
            r_subset_p = _is_multiset_subset(r_frags, p_frags)
            p_subset_r = _is_multiset_subset(p_frags, r_frags)
            if r_subset_p or p_subset_r:
                final_reactant_match = reactant_matches[i]
                final_product_match = product_matches[j]
                break
        if final_reactant_match is not None:
            break

    # No acceptable pairing found
    if final_reactant_match is None or final_product_match is None:
        return {
            "rxn": rxn,
            "mapped_atoms": mapped_atoms,
            "unmapped_matched_atoms": unmapped_matched_atoms,
            "unmatched_atoms": unmatched_atoms,
        }

    # Classify atoms on reactant
    for i, rxn_idx in enumerate(final_reactant_match):
        op_atom = operator_reactant.GetAtomWithIdx(i)
        map_num = op_atom.GetAtomMapNum()
        if map_num > 0:
            mapped_atoms[0].append((rxn_idx, map_num))
        else:
            unmapped_matched_atoms[0].append(rxn_idx)

    # Classify atoms on product
    for i, rxn_idx in enumerate(final_product_match):
        op_atom = operator_product.GetAtomWithIdx(i)
        map_num = op_atom.GetAtomMapNum()
        if map_num > 0:
            mapped_atoms[1].append((rxn_idx, map_num))
        else:
            unmapped_matched_atoms[1].append(rxn_idx)

    # Unmatched atoms are those not in mapped or unmapped_matched
    all_r_indices = set(range(rxn_reactant.GetNumAtoms()))
    all_p_indices = set(range(rxn_product.GetNumAtoms()))

    matched_r = set(idx for idx, _ in mapped_atoms[0]) | set(unmapped_matched_atoms[0])
    matched_p = set(idx for idx, _ in mapped_atoms[1]) | set(unmapped_matched_atoms[1])

    unmatched_atoms = (
        sorted(all_r_indices - matched_r),
        sorted(all_p_indices - matched_p),
    )

    return {
        "rxn": rxn,
        "mapped_atoms": mapped_atoms,
        "unmapped_matched_atoms": unmapped_matched_atoms,
        "unmatched_atoms": unmatched_atoms,
    }


def _has_nontrivial_match(result: Dict) -> bool:
    """
    Return True if the mechanistic labeling contains any matched atoms
    (mapped or unmapped). Used to decide whether an operator matches
    a concrete reaction.
    """
    mapped_r, mapped_p = result["mapped_atoms"]
    unmapped_r, unmapped_p = result["unmapped_matched_atoms"]
    return bool(mapped_r or mapped_p or unmapped_r or unmapped_p)


def _reduced_fragment_smiles(mol: Chem.Mol, remove_indices) -> List[str]:
    """
    Remove atoms in remove_indices, then return canonical SMILES for each
    disconnected fragment of the remainder.
    """
    rw_mol = Chem.RWMol()
    idx_map = {}

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in remove_indices:
            new_idx = rw_mol.AddAtom(Chem.Atom(atom.GetAtomicNum()))
            idx_map[idx] = new_idx

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in idx_map and a2 in idx_map:
            rw_mol.AddBond(idx_map[a1], idx_map[a2], bond.GetBondType())

    reduced = rw_mol.GetMol()
    if reduced.GetNumAtoms() == 0:
        return []

    frags = Chem.GetMolFrags(reduced, asMols=True)
    return [Chem.MolToSmiles(f, canonical=True) for f in frags]


def _is_multiset_subset(a: List[str], b: List[str]) -> bool:
    """
    Return True if multiset a is a subset of multiset b, using fragment SMILES.
    """
    ca = Counter(a)
    cb = Counter(b)
    return all(ca[k] <= cb[k] for k in ca)


def _indexed_smiles(mol: Chem.Mol) -> str:
    """
    Return a SMILES string where each atom's map number is its RDKit atom index.
    This makes it easy to interpret integer indices in the labeling output.
    """
    tmp = Chem.Mol(mol)  # shallow copy
    for atom in tmp.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return Chem.MolToSmiles(tmp, canonical=True)


def _format_atom_list(mol: Chem.Mol, indices: List[int]) -> str:
    if not indices:
        return "none"
    return ", ".join(f"{idx} ({mol.GetAtomWithIdx(idx).GetSymbol()})" for idx in indices)


def _pretty_print_mechanistic_result(result: Dict) -> None:
    rxn = result["rxn"]
    rxn_reactant = rxn.GetReactants()[0]
    rxn_product = rxn.GetProducts()[0]

    mapped_r, mapped_p = result["mapped_atoms"]
    unmapped_r, unmapped_p = result["unmapped_matched_atoms"]
    unmatched_r, unmatched_p = result["unmatched_atoms"]

    print("=== Mechanistic labeling summary ===")
    print("Substrate (SMILES):", Chem.MolToSmiles(rxn_reactant))
    print("Product   (SMILES):", Chem.MolToSmiles(rxn_product))
    print("Substrate (indexed SMILES; atom_idx as map):", _indexed_smiles(rxn_reactant))
    print("Product   (indexed SMILES; atom_idx as map):", _indexed_smiles(rxn_product))
    print()

    # Organize mapped atoms by map number
    r_by_map = {}
    for idx, map_num in mapped_r:
        r_by_map.setdefault(map_num, []).append(idx)
    p_by_map = {}
    for idx, map_num in mapped_p:
        p_by_map.setdefault(map_num, []).append(idx)

    all_map_nums = sorted(set(r_by_map.keys()) | set(p_by_map.keys()))

    print("Mapped atoms (by map number):")
    if not all_map_nums:
        print("  none")
    else:
        for map_num in all_map_nums:
            r_idxs = r_by_map.get(map_num, [])
            p_idxs = p_by_map.get(map_num, [])
            r_descr = _format_atom_list(rxn_reactant, r_idxs)
            p_descr = _format_atom_list(rxn_product, p_idxs)
            print(f"  map {map_num}: reactant {r_descr} -> product {p_descr}")
    print()

    print("Unmapped but matched atoms:")
    print("  Reactant:", _format_atom_list(rxn_reactant, unmapped_r))
    print("  Product: ", _format_atom_list(rxn_product, unmapped_p))
    print()

    print("Unmatched atoms:")
    print("  Reactant:", _format_atom_list(rxn_reactant, unmatched_r))
    print("  Product: ", _format_atom_list(rxn_product, unmatched_p))
    print()


def find_mechanistic_match_in_dataset(
    reaction_smiles: str,
    abstraction: str,
    data_dir: Optional[Path] = None,
) -> Optional[Dict[str, Any]]:
    """
    Scan an EVODEX mechanistic operator dataset (Bm, Cm, Dm, Em) from top to
    bottom and return the first operator whose mechanistic labeling matches
    the given reaction.

    Parameters
    ----------
    reaction_smiles : str
        Concrete reaction SMILES to be labeled (reactant>>product).
        Must be a single substrate to single product reaction.
    abstraction : str
        Mechanistic abstraction level. Expected values are "Bm", "Cm", "Dm",
        or "Em". Case insensitive.
    data_dir : pathlib.Path, optional
        Directory containing the EVODEX CSV files. Defaults to the "data"
        subdirectory next to this module.

    Returns
    -------
    dict or None
        On success, returns a dict with keys:

          - "dataset_id": EVODEX operator id (for example "EVODEX.2-Em1").
          - "operator_smirks": operator SMIRKS from the CSV.
          - "sources": sources field from the CSV (if present).
          - "hash": hash field from the CSV (if present).
          - "labeling": full mechanistic labeler result dict.

        If no operator matches the reaction, returns None.

    Raises
    ------
    ValueError
        If the reaction SMILES contains multiple molecules (indicated by '.').
    """
    abstraction_clean = abstraction.strip().lower()

    file_map = {
        "bm": "EVODEX-Bm.csv",
        "cm": "EVODEX-Cm.csv",
        "dm": "EVODEX-Dm.csv",
        "em": "EVODEX-Em.csv",
    }

    if abstraction_clean not in file_map:
        raise ValueError(
            f"Unsupported abstraction level {abstraction!r}. "
            f"Expected one of {sorted(file_map.keys())} (case insensitive)."
        )

    if data_dir is None:
        data_dir = Path(__file__).resolve().parent / "data"

    csv_path = data_dir / file_map[abstraction_clean]

    if not csv_path.is_file():
        raise FileNotFoundError(f"EVODEX dataset not found at {csv_path}")

    # Prepare the concrete reaction once (will raise ValueError if multi-molecule)
    rxn = prepare_reaction(reaction_smiles)

    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smirks = row.get("smirks")
            if not smirks:
                continue
            
            # Skip multi-molecule operators
            if '.' in smirks:
                continue

            operator = prepare_operator(smirks)
            result = mechanistic_label_reaction(rxn, operator)

            if _has_nontrivial_match(result):
                return {
                    "dataset_id": row.get("id"),
                    "operator_smirks": smirks,
                    "sources": row.get("sources"),
                    "hash": row.get("hash"),
                    "labeling": result,
                }

    # No match found
    return None


if __name__ == "__main__":
    # Verbose example: kinase like complete operator with explicit hydrogens.
    # This is meant for interactive inspection, not for unit tests.
    operator_smirks = "[C:1][O:2][H]>>[C:1][O:2]P(=O)(O[H])O[H]"
    reaction_smiles = "CCCO>>CCCOP(=O)(O)O"

    print("Operator SMIRKS:", operator_smirks)
    print("Reaction SMILES:", reaction_smiles)

    # Prepare reaction with explicit hydrogens, matching the complete operator style.
    rxn = prepare_reaction(reaction_smiles)
    rxn_reactant = rxn.GetReactants()[0]
    rxn_product = rxn.GetProducts()[0]

    print("\nPrepared Reaction (SMILES with explicit Hs):")
    print(rdChemReactions.ReactionToSmiles(rxn))

    operator = prepare_operator(operator_smirks)
    op_r = operator.GetReactants()[0]
    op_p = operator.GetProducts()[0]

    print("\nOperator templates as RDKit SMARTS:")
    print("  Reactant:", Chem.MolToSmarts(op_r))
    print("  Product :", Chem.MolToSmarts(op_p))

    # Finally, run the mechanistic labeler and show the human readable summary.
    result = mechanistic_label_reaction(rxn, operator)

    print("\n=== Mechanistic labeler result ===")
    _pretty_print_mechanistic_result(result)

    # Example of scanning the Em dataset for a mechanistic match
    # (requires EVODEX-Em.csv in evodex/data)
    try:
        match = find_mechanistic_match_in_dataset(reaction_smiles, "Em")
        if match is not None:
            print("\nFirst Em dataset match:")
            print("  id:", match["dataset_id"])
            print("  operator:", match["operator_smirks"])
        else:
            print("\nNo Em dataset match found for this reaction.")
    except FileNotFoundError as e:
        print("\nSkipping dataset scan:", e)