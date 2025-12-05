from typing import Dict, Tuple, List
from collections import Counter

from rdkit import Chem
from rdkit.Chem import rdChemReactions, AllChem


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

    This helper:
      - Splits the reaction SMILES into reactant and product parts.
      - Parses individual molecules on each side.
      - Adds explicit hydrogens to all atoms.
      - Constructs a ChemicalReaction with one or more reactant and product templates.

    Note
    ----
    Downstream labelers such as mechanistic_label_reaction assume a
    single-substrate and single-product viewpoint and operate on
    rxn.GetReactants()[0] and rxn.GetProducts()[0]. This helper does not
    enforce that there is exactly one reactant or product template, so callers
    that rely on a 1-to-1 substrate to product mapping must handle that
    assumption explicitly.

    Parameters
    ----------
    reaction_smiles : str
        A reaction SMILES string of the form 'reactants>>products'.

    Returns
    -------
    rdkit.Chem.rdChemReactions.ChemicalReaction
        An RDKit ChemicalReaction object with hydrogens added to all
        reactant and product molecules.
    """
    reactant_smiles, product_smiles = reaction_smiles.split(">>")
    reactant_mols = [
        Chem.AddHs(Chem.MolFromSmiles(smi))
        for smi in reactant_smiles.split(".")
        if smi
    ]
    product_mols = [
        Chem.AddHs(Chem.MolFromSmiles(smi))
        for smi in product_smiles.split(".")
        if smi
    ]

    rxn = AllChem.ChemicalReaction()
    for mol in reactant_mols:
        rxn.AddReactantTemplate(mol)
    for mol in product_mols:
        rxn.AddProductTemplate(mol)

    return rxn


def mechanistic_label_reaction(rxn, operator) -> Dict:
    """
    Perform a relaxed, mechanistically oriented alignment between a concrete reaction
    and a reaction operator.

    This function conceptually treats the reaction as a single-substrate and
    single-product transformation and uses the first reactant and first product
    templates in each ChemicalReaction. It identifies atom-level
    correspondences between the operator and the reaction and classifies
    reaction atoms by their parity in the transformation.

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
