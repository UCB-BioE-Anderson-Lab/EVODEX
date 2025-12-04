from typing import Set, Tuple

import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from evodex.mechanistic_labeler import (
    prepare_reaction,
    prepare_operator,
    mechanistic_label_reaction,
)


def _assert_empty_label(result: dict) -> None:
    """Helper to assert the canonical 'no match' structure."""
    assert result["mapped_atoms"] == ([], [])
    assert result["unmapped_matched_atoms"] == ([], [])
    assert result["unmatched_atoms"] == ([], [])


def _collect_atom_index_partition(result: dict) -> Tuple[Set[int], Set[int], Set[int]]:
    """
    Return sets of indices (reactant only) for mapped, unmapped-matched, and unmatched.
    This is used only for partition invariants, not for checking specific indices.
    """
    mapped_r = {idx for idx, _ in result["mapped_atoms"][0]}
    unmapped_r = set(result["unmapped_matched_atoms"][0])
    unmatched_r = set(result["unmatched_atoms"][0])
    return mapped_r, unmapped_r, unmatched_r


def test_prepare_reaction_single_templates() -> None:
    """prepare_reaction should produce a 1→1 reaction with explicit hydrogens."""
    reaction_smiles = "OCC(C)CCCCO>>OCC(C)CCCC=O"
    rxn = prepare_reaction(reaction_smiles)

    assert rxn.GetNumReactantTemplates() == 1
    assert rxn.GetNumProductTemplates() == 1

    # Sanity: adding Hs should increase atom counts relative to bare SMILES
    reactant = rxn.GetReactants()[0]
    bare = Chem.MolFromSmiles("OCC(C)CCCCO")
    assert reactant.GetNumAtoms() > bare.GetNumAtoms()


def test_mechanistic_label_reaction_returns_same_rxn() -> None:
    """mechanistic_label_reaction must return the same rxn object it was given."""
    reaction_smiles = "OCC(C)CCCCO>>OCC(C)CCCC=O"
    rxn = prepare_reaction(reaction_smiles)
    operator_smarts = "[C:1]([H])[O:2]([H])>>[C:1]=[O:2]"
    operator = rdChemReactions.ReactionFromSmarts(operator_smarts)

    result = mechanistic_label_reaction(rxn, operator)

    assert result["rxn"] is rxn


def test_mechanistic_label_redox_equal_fragments() -> None:
    """
    Primary alcohol oxidation example (same as original label_reaction plumbing).

    Here the non-transformed context is identical up to explicit H decoration,
    so the mechanistic labeler should find a non-empty mapping whose map-number
    sets exactly match the operator.
    """
    reaction_smiles = "OCC(C)CCCCO>>OCC(C)CCCC=O"
    rxn = prepare_reaction(reaction_smiles)
    operator_smarts = "[C:1]([H])[O:2]([H])>>[C:1]=[O:2]"
    operator = rdChemReactions.ReactionFromSmarts(operator_smarts)

    result = mechanistic_label_reaction(rxn, operator)

    mapped_r, mapped_p = result["mapped_atoms"]
    assert mapped_r
    assert mapped_p

    # Operator maps C:1 and O:2 on both sides.
    maps_r = {m for _, m in mapped_r}
    maps_p = {m for _, m in mapped_p}
    assert maps_r == {1, 2}
    assert maps_p == {1, 2}


def test_mechanistic_label_esterase_reactant_superset() -> None:
    """
    Esterase-like case: extra context on the reactant side (2→1 fragment scenario).

    Operator only describes the leaving alcohol fragment; the acyl fragment remains
    unmatched context on the reactant.
    """
    operator_smirks = "[C:1](=O)[O:2][C:3]>>[H][O:2][C:3]"
    substrate_smiles = "CCOC(=O)CCC"
    product_smiles = "CCO"

    reaction_smiles = f"{substrate_smiles}>>{product_smiles}"
    rxn = prepare_reaction(reaction_smiles)
    operator = prepare_operator(operator_smirks)

    result = mechanistic_label_reaction(rxn, operator)

    mapped_r, mapped_p = result["mapped_atoms"]
    assert mapped_r
    assert mapped_p

    # Operator definition:
    #   - Reactant maps C:1, O:2, C:3
    #   - Product maps O:2, C:3
    maps_r = {m for _, m in mapped_r}
    maps_p = {m for _, m in mapped_p}
    assert maps_r == {1, 2, 3}
    assert maps_p == {2, 3}

    unmatched_r, unmatched_p = result["unmatched_atoms"]
    # Extra acyl context: reactant must have strictly more unmatched atoms.
    assert len(unmatched_r) > len(unmatched_p)


def test_mechanistic_label_kinase_product_superset_complete_operator() -> None:
    """
    Kinase-like toy case using a complete operator with explicit hydrogens.

    The operator describes phosphorylation of an alcohol including the O–H
    protons; the phosphate context appears as matched but (mostly) unmapped
    atoms on the product side.
    """
    operator_smirks = "[C:1][O:2][H]>>[C:1][O:2]P(=O)(O[H])O[H]"
    substrate_smiles = "CCO"
    product_smiles = "CCOP(=O)(O)O"

    reaction_smiles = f"{substrate_smiles}>>{product_smiles}"
    rxn = prepare_reaction(reaction_smiles)
    operator = prepare_operator(operator_smirks)

    result = mechanistic_label_reaction(rxn, operator)

    mapped_r, mapped_p = result["mapped_atoms"]

    assert mapped_r
    assert mapped_p

    # Operator maps C:1 and O:2 on both sides.
    maps_r = {m for _, m in mapped_r}
    maps_p = {m for _, m in mapped_p}
    assert maps_r == {1, 2}
    assert maps_p == {1, 2}

    unmapped_r, unmapped_p = result["unmapped_matched_atoms"]
    # Complete operator should introduce unmapped matched atoms
    # (e.g. O–H protons / phosphate context) at least on the product side.
    assert unmapped_p


def test_mechanistic_label_kinase_product_superset_stripped_operator() -> None:
    """
    Kinase-like toy case with a stripped operator (no explicit O–H hydrogens).

    Here the phosphate group is treated as extra context on the product side,
    so there should be strictly more unmatched atoms on the product.
    """
    operator_smirks = "[C:1][O:2][H]>>[C:1][O:2][P]"
    substrate_smiles = "CCO"
    product_smiles = "CCOP(=O)(O)O"

    reaction_smiles = f"{substrate_smiles}>>{product_smiles}"
    rxn = prepare_reaction(reaction_smiles)
    operator = prepare_operator(operator_smirks)

    result = mechanistic_label_reaction(rxn, operator)
    print(result)

    mapped_r, mapped_p = result["mapped_atoms"]
    assert mapped_r
    assert mapped_p

    # Same mapping set on both sides (C:1 and O:2).
    maps_r = {m for _, m in mapped_r}
    maps_p = {m for _, m in mapped_p}
    assert maps_r == {1, 2}
    assert maps_p == {1, 2}

    unmatched_r, unmatched_p = result["unmatched_atoms"]
    # Extra phosphate context on the product side.
    assert len(unmatched_p) > len(unmatched_r)


def test_mechanistic_label_no_reactant_match_returns_empty() -> None:
    """If the operator does not match the substrate at all, the labeler returns empty labels."""
    operator_smirks = "[N:1]>>[N:1]"
    substrate_smiles = "CCC"
    product_smiles = "CCC"

    reaction_smiles = f"{substrate_smiles}>>{product_smiles}"
    rxn = prepare_reaction(reaction_smiles)
    operator = prepare_operator(operator_smirks)

    result = mechanistic_label_reaction(rxn, operator)

    _assert_empty_label(result)


def test_mechanistic_label_no_product_match_returns_empty() -> None:
    """If the operator does not match the product at all, the labeler returns empty labels."""
    operator_smirks = "[C:1]([H])[O:2]([H])>>[C:1]=[O:2]"
    substrate_smiles = "CCO"
    product_smiles = "CCC"  # no carbonyl

    reaction_smiles = f"{substrate_smiles}>>{product_smiles}"
    rxn = prepare_reaction(reaction_smiles)
    operator = prepare_operator(operator_smirks)

    result = mechanistic_label_reaction(rxn, operator)

    _assert_empty_label(result)


def test_mechanistic_label_no_pairing_returns_empty() -> None:
    """
    Both sides match individually, but no reactant/product pairing satisfies the fragment subset rule.

    Here the operator matches a CO fragment; substrates differ in the leftover fragment
    (Cl vs Br), so the fragment multisets are incomparable.
    """
    operator_smirks = "[C:1][O:2]>>[C:1][O:2]"
    substrate_smiles = "COCl"
    product_smiles = "COBr"

    reaction_smiles = f"{substrate_smiles}>>{product_smiles}"
    rxn = prepare_reaction(reaction_smiles)
    operator = prepare_operator(operator_smirks)

    result = mechanistic_label_reaction(rxn, operator)

    _assert_empty_label(result)


def test_mechanistic_label_atom_partition_is_disjoint_and_complete() -> None:
    """
    For a simple successful case, the union of mapped, unmapped-matched, and unmatched
    atom indices on the reactant side must cover all atoms and be pairwise disjoint.
    """
    reaction_smiles = "OCC(C)CCCCO>>OCC(C)CCCC=O"
    rxn = prepare_reaction(reaction_smiles)
    operator_smarts = "[C:1]([H])[O:2]([H])>>[C:1]=[O:2]"
    operator = rdChemReactions.ReactionFromSmarts(operator_smarts)

    result = mechanistic_label_reaction(rxn, operator)

    mapped_r, unmapped_r, unmatched_r = _collect_atom_index_partition(result)

    reactant = rxn.GetReactants()[0]
    all_indices = set(range(reactant.GetNumAtoms()))

    # Pairwise disjoint.
    assert mapped_r.isdisjoint(unmapped_r)
    assert mapped_r.isdisjoint(unmatched_r)
    assert unmapped_r.isdisjoint(unmatched_r)

    # Complete cover.
    assert mapped_r | unmapped_r | unmatched_r == all_indices
