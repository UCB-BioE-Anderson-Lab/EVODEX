"""Utilities for mapping explicit hydrogen atoms in reaction SMILES.

This module takes a reaction SMILES string, parses reactants and products,
adds explicit hydrogens and assigns atom maps to many of the hydrogens
in a way that is consistent between reactants and products.

It begins with a heavy-atom-mapped reaction SMILES without explicit hydrgens
and thus no hydrogen atom maps.  It identifies the hydrogen maps by the
consistency of atom mapping on the attached heavy atom.  If that heavy atom
has a hydrogen on both the substrate and product, the hydrogens get the next
available atom map.  If the heavy atom only has the hydrogen on one side of
the reaction, it remains unmapped.  Thus, even in a full balanced reaction
there will be hydrogens that remain unmapped. When constructing operators from
these, the lack of parity of these hydrogens implies that it is a reaction
center and will be incorporated into the final operator at all levels of
evodex abstractions. Thus, it does not ultimately matter that they are
unmapped.

High level algorithm
--------------------
1. Split the reaction SMILES over '>>' into reactant and product strings.
2. Split each side over '.' to get individual molecule SMILES.
3. Convert each molecule to an RDKit Mol and call Chem.AddHs.
4. Build a ChemicalReaction from those hydrogenated molecules.
5. Validate that *heavy* atoms (non‑hydrogen) already have unique,
   positive atom‑map numbers and that the sets of heavy‑atom maps
   are identical between reactants and products.
6. Walk over all explicit hydrogens on the reactant side, assigning
   new map numbers and recording a mapping keyed by the heavy atom
   they are attached to.
7. Walk over explicit hydrogens on the product side and, whenever
   possible, re‑use the hydrogen map numbers that were created for
   the corresponding heavy‑atom maps on the reactant side.
8. Any remaining hydrogens that could not be consistently matched
   are left unmapped (map number 0).

The main public entry point is :func:`map_hydrogens_in_reaction` which
returns a new reaction SMILES string with explicit hydrogens and the
assigned atom maps.
"""

from __future__ import annotations

from typing import Iterable, List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import ReactionToImage
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def _indexed_smiles(mol: Chem.Mol) -> str:
    """Return a SMILES string where each atom's map number is its RDKit index.

    This is primarily a debugging helper – it makes it easy to interpret
    integer indices in labeling output, visualisation tools, or tests.
    """
    tmp = Chem.Mol(mol)  # shallow copy
    for atom in tmp.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return Chem.MolToSmiles(tmp, canonical=True)


def _split_reaction_smiles(rxn_smiles: str) -> Tuple[List[str], List[str]]:
    """Split a reaction SMILES into lists of reactant and product SMILES.

    Parameters
    ----------
    rxn_smiles:
        Reaction SMILES of the form "reactants>>products". Each side can
        contain multiple molecules separated by '.'.
    """
    if ">>" not in rxn_smiles:
        raise ValueError(f"Reaction SMILES must contain '>>': {rxn_smiles}")

    reactant_part, product_part = rxn_smiles.split(">>", 1)

    def _split_side(side: str) -> List[str]:
        side = side.strip()
        if not side:
            return []
        return [s for s in (p.strip() for p in side.split(".")) if s]

    reactants = _split_side(reactant_part)
    products = _split_side(product_part)
    return reactants, products


def _mols_from_smiles(smiles_list: Iterable[str]) -> List[Chem.Mol]:
    """Convert an iterable of SMILES strings to RDKit molecules.

    Raises a ValueError if any SMILES cannot be parsed.
    """
    mols: List[Chem.Mol] = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Could not parse SMILES: {smi}")
        mols.append(mol)
    return mols


def _add_hydrogens(mols: Iterable[Chem.Mol]) -> List[Chem.Mol]:
    """Return copies of the input molecules with explicit hydrogens added."""
    return [Chem.AddHs(mol) for mol in mols]


def _build_reaction_with_hydrogens(rxn_smiles: str) -> AllChem.ChemicalReaction:
    """Parse a reaction SMILES and return a ChemicalReaction with explicit Hs.

    The heavy atoms retain whatever atom‑map numbers are present in the
    input SMILES. Explicit hydrogens are added with map number 0.
    """
    reactant_smis, product_smis = _split_reaction_smiles(rxn_smiles)

    reactant_mols = _add_hydrogens(_mols_from_smiles(reactant_smis))
    product_mols = _add_hydrogens(_mols_from_smiles(product_smis))

    rxn = AllChem.ChemicalReaction()
    for mol in reactant_mols:
        rxn.AddReactantTemplate(mol)
    for mol in product_mols:
        rxn.AddProductTemplate(mol)

    return rxn


def _validate_heavy_atom_maps(mol: Chem.Mol, smiles: str) -> set[int]:
    """Validate atom‑map numbers for heavy atoms in a molecule.

    Hydrogens (atomic number 1) are ignored here because their maps are
    assigned separately in :func:`map_hydrogens_in_reaction`.

    Returns
    -------
    set[int]
        The set of heavy‑atom map numbers present in the molecule.
    """
    atom_map_set: set[int] = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:  # skip hydrogens
            continue
        atom_map = atom.GetAtomMapNum()
        if atom_map <= 0:
            raise ValueError(
                f"Heavy atom without valid map number: {atom.GetSymbol()} in {smiles}"
            )
        if atom_map in atom_map_set:
            raise ValueError(
                f"Duplicate atom map number {atom_map} in heavy atoms for {smiles}"
            )
        atom_map_set.add(atom_map)
    return atom_map_set


def _collect_heavy_atom_maps(rxn: AllChem.ChemicalReaction, smiles: str) -> Tuple[set[int], set[int]]:
    """Return sets of heavy‑atom map numbers for reactants and products."""
    reactant_maps: set[int] = set()
    product_maps: set[int] = set()

    for reactant in rxn.GetReactants():
        reactant_maps.update(_validate_heavy_atom_maps(reactant, smiles))

    for product in rxn.GetProducts():
        product_maps.update(_validate_heavy_atom_maps(product, smiles))

    return reactant_maps, product_maps


def map_hydrogens_in_reaction(rxn_smiles: str) -> str:
    """Assign atom maps to many explicit hydrogens in a reaction SMILES.

    Parameters
    ----------
    rxn_smiles:
        Reaction SMILES string. Heavy atoms are expected to already carry
        valid, unique atom‑map numbers that are consistent between
        reactants and products. Hydrogens may or may not be explicit; any
        that are present will be made explicit and (where possible) given
        map numbers that are matched between reactants and products.

    Returns
    -------
    str
        A reaction SMILES string with explicit hydrogens and updated atom
        maps.
    """
    # 1–4: parse, split, build a ChemicalReaction with explicit hydrogens
    try:
        reaction = _build_reaction_with_hydrogens(rxn_smiles)
    except Exception as exc:  # keep the original context in the message
        raise ValueError(f"Invalid reaction SMILES: {rxn_smiles}") from exc

    # 5: validate heavy‑atom maps and enforce consistency
    reactant_maps, product_maps = _collect_heavy_atom_maps(reaction, rxn_smiles)

    if reactant_maps != product_maps:
        raise ValueError(
            "Mismatch between reactant and product heavy‑atom maps in "
            f"{rxn_smiles}"
        )

    if not reactant_maps:
        # Degenerate case – nothing is mapped; return a hydrogenated version
        # without attempting to assign hydrogen maps.
        reactant_smiles = [
            Chem.MolToSmiles(mol, isomericSmiles=True)
            for mol in reaction.GetReactants()
        ]
        product_smiles = [
            Chem.MolToSmiles(mol, isomericSmiles=True)
            for mol in reaction.GetProducts()
        ]
        return ">>".join([".".join(reactant_smiles), ".".join(product_smiles)])

    next_atom_map = max(reactant_maps.union(product_maps)) + 1

    # For each heavy‑atom map number, retain the list of hydrogen map
    # numbers that were assigned on the reactant side. These are reused
    # for hydrogens attached to the same heavy‑atom map in the products.
    heavy_to_h_map: dict[int, List[int]] = {}

    def _assign_reactant_h_maps(mol: Chem.Mol) -> None:
        nonlocal next_atom_map
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            neighbors = list(atom.GetNeighbors())
            if not neighbors:
                continue
            # Use the first neighbor – for normal covalent hydrogens this
            # is unambiguous.
            neighbor = neighbors[0]
            neighbor_map = neighbor.GetAtomMapNum()
            if neighbor_map <= 0:
                # No meaningful heavy‑atom map to attach to; leave unmapped
                continue

            if neighbor_map not in heavy_to_h_map:
                heavy_to_h_map[neighbor_map] = []

            atom.SetAtomMapNum(next_atom_map)
            heavy_to_h_map[neighbor_map].append(next_atom_map)
            next_atom_map += 1

    def _assign_product_h_maps(mol: Chem.Mol) -> None:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            neighbors = list(atom.GetNeighbors())
            if not neighbors:
                continue
            neighbor = neighbors[0]
            neighbor_map = neighbor.GetAtomMapNum()
            if neighbor_map <= 0:
                continue

            queue = heavy_to_h_map.get(neighbor_map)
            if queue:
                # Reuse an existing hydrogen map from the reactant side
                atom.SetAtomMapNum(queue.pop(0))
            else:
                # No corresponding reactant hydrogen; leave this unmapped
                atom.SetAtomMapNum(0)

    # 6: assign hydrogen maps on the reactant side
    for reactant in reaction.GetReactants():
        _assign_reactant_h_maps(reactant)

    # 7: assign hydrogen maps on the product side, reusing where possible
    for product in reaction.GetProducts():
        _assign_product_h_maps(product)

    # 8: any hydrogen on the reactant side whose map number is still present
    # in the queues was never matched to a product hydrogen. To avoid
    # dangling, one‑sided maps, clear their map number back to zero.
    remaining_maps = {
        h_map
        for queue in heavy_to_h_map.values()
        for h_map in queue
    }

    if remaining_maps:
        remaining_maps = set(remaining_maps)
        for reactant in reaction.GetReactants():
            for atom in reactant.GetAtoms():
                if atom.GetAtomicNum() == 1 and atom.GetAtomMapNum() in remaining_maps:
                    atom.SetAtomMapNum(0)

    # Convert back to a reaction SMILES string. We use SMARTS writers to
    # preserve atom maps as they appear in the reaction object.
    try:
        reactant_smiles = [
            Chem.MolToSmarts(mol, isomericSmiles=True)
            for mol in reaction.GetReactants()
        ]
        product_smiles = [
            Chem.MolToSmarts(mol, isomericSmiles=True)
            for mol in reaction.GetProducts()
        ]
        modified_smiles = ">>".join([
            ".".join(reactant_smiles),
            ".".join(product_smiles),
        ])
    except Exception as exc:
        raise ValueError(
            f"Error converting modified reaction to SMIRKS for {rxn_smiles}"
        ) from exc

    return modified_smiles


__all__ = [
    "_indexed_smiles",
    "map_hydrogens_in_reaction",
]

if __name__ == "__main__":
    example_rxn = (
        # "[CH3:1][S+:2]([CH2:3][CH2:4][C@H:5]([NH2:6])[C:7](=[O:8])[OH:9])[CH2:10][C@H:11]1[O:12][C@@H:13]([n:14]2[cH:15][n:16][c:17]3[c:18]([NH2:19])[n:20][cH:21][n:22][c:23]23)[C@H:24]([OH:25])[C@@H:26]1[OH:27].[O:28]=[C:29]1[CH2:30][C@@H:31]([c:32]2[cH:33][cH:34][c:35]([OH:36])[c:37]([OH:38])[cH:39]2)[O:40][c:41]2[cH:42][c:43]([OH:44])[cH:45][c:46]([OH:47])[c:48]21>>[CH3:1][O:38][c:37]1[c:35]([OH:36])[cH:34][cH:33][c:32]([C@@H:31]2[CH2:30][C:29](=[O:28])[c:48]3[c:41]([cH:42][c:43]([OH:44])[cH:45][c:46]3[OH:47])[O:40]2)[cH:39]1.[H+].[S:2]([CH2:3][CH2:4][C@H:5]([NH2:6])[C:7](=[O:8])[OH:9])[CH2:10][C@H:11]1[O:12][C@@H:13]([n:14]2[cH:15][n:16][c:17]3[c:18]([NH2:19])[n:20][cH:21][n:22][c:23]23)[C@H:24]([OH:25])[C@@H:26]1[OH:27]",
        "[O:1]=[C:2]([OH:3])[CH2:4][C:5](=[O:6])[C:7](=[O:8])[OH:9]>>[CH3:4][C:5](=[O:6])[C:7](=[O:8])[OH:9].[O:1]=[C:2]=[O:3]"
    
    )
    print("Input reaction SMILES:")
    print(example_rxn)
    print()
    print("Reaction with explicit hydrogen maps:")
    mapped_rxn = map_hydrogens_in_reaction(example_rxn)
    print(mapped_rxn)

    # Read the mapped reaction back in and generate an image
    try:
        rxn_obj = AllChem.ReactionFromSmarts(mapped_rxn, useSmiles=True)
        img = ReactionToImage(rxn_obj)
        img.save("hydrogen_mapped_example.png")
        print()
        print("Saved reaction image to hydrogen_mapped_example.png")
    except Exception as exc:
        print()
        print(f"Could not generate reaction image: {exc}")
