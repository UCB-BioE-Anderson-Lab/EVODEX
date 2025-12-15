#!/usr/bin/env python3

import csv
import os
from typing import Set, Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import inchi
from rdkit import RDLogger

from evodex.utils import get_molecule_hash_from_mol

RDLogger.DisableLog("rdApp.*")

# Global cache for native metabolites (hashed)
_native_metabolites: Optional[Set[str]] = None


def _load_native_metabolites() -> Set[str]:
    """
    Load ubiquitous metabolites from a TSV and return a set of canonical hashes.
    Hashing is done directly from RDKit Mol via get_molecule_hash_from_mol.
    """
    global _native_metabolites
    if _native_metabolites is not None:
        return _native_metabolites

    script_dir = os.path.dirname(__file__)
    ubiquitous_metabolites_file = os.path.join(
        script_dir, "data", "ubiquitous_metabolites.txt"
    )
    backup_file = os.path.join(
        script_dir, "data", "ubiquitous_metabolites_backup.txt"
    )

    native_hashes: Set[str] = set()

    try:
        with open(ubiquitous_metabolites_file, "r", newline="") as file:
            reader = csv.DictReader(file, delimiter="\t")

            with open(backup_file, "w", newline="") as backup:
                backup_writer = csv.writer(backup, delimiter="\t")
                backup_writer.writerow(["name", "original_inchi", "custom_hash"])

                for row in reader:
                    name = (row.get("name") or "").strip().strip('"')
                    inchi_str = (row.get("inchi") or "").strip().strip('"')

                    if not inchi_str:
                        continue

                    mol = inchi.MolFromInchi(inchi_str)
                    if mol is None:
                        continue

                    try:
                        # Hash directly from Mol (no SMILES round-trip)
                        custom_hash = get_molecule_hash_from_mol(mol)
                        native_hashes.add(custom_hash)
                        backup_writer.writerow([name, inchi_str, custom_hash])
                    except Exception as e:
                        print(f"Failed to process molecule: {name}, Error: {e}")

    except Exception as e:
        raise RuntimeError(f"Failed to load ubiquitous metabolites: {e}")

    _native_metabolites = native_hashes
    return _native_metabolites


def _clean_up_atom_maps(rxn: AllChem.ChemicalReaction) -> None:
    """
    Remove atom maps in products that do not originate from any reactant atom map.
    Also removes dangling reactant atom maps that never appear in products.
    """
    try:
        substrate_atom_maps = set()

        for mol in rxn.GetReactants():
            for atom in mol.GetAtoms():
                m = atom.GetAtomMapNum()
                if m > 0:
                    substrate_atom_maps.add(m)

        for mol in rxn.GetProducts():
            for atom in mol.GetAtoms():
                m = atom.GetAtomMapNum()
                if m > 0:
                    if m not in substrate_atom_maps:
                        atom.SetAtomMapNum(0)
                    else:
                        substrate_atom_maps.remove(m)

        for mol in rxn.GetReactants():
            for atom in mol.GetAtoms():
                m = atom.GetAtomMapNum()
                if m in substrate_atom_maps:
                    atom.SetAtomMapNum(0)

    except Exception as e:
        raise RuntimeError(f"Failed to clean up atom maps: {e}")


def remove_cofactors(rxn_smiles: str) -> str:
    """
    Given a reaction SMILES/SMIRKS-like string (parsed with useSmiles=True),
    remove any reactant/product molecules whose hash is in the ubiquitous metabolite set.
    Returns a SMIRKS-like string composed from RDKit MolToSmarts on the filtered templates.
    """
    try:
        native_metabolites = _load_native_metabolites()

        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            raise ValueError(f"Invalid reaction string: {rxn_smiles}")

        non_cofactor_reactants = []
        non_cofactor_products = []

        for mol in rxn.GetReactants():
            try:
                custom_hash = get_molecule_hash_from_mol(mol)
                if custom_hash not in native_metabolites:
                    non_cofactor_reactants.append(
                        Chem.MolToSmiles(mol, isomericSmiles=True)
                    )
            except Exception as e:
                raise RuntimeError(f"Failed to process reactant: {e}")

        for mol in rxn.GetProducts():
            try:
                custom_hash = get_molecule_hash_from_mol(mol)
                if custom_hash not in native_metabolites:
                    non_cofactor_products.append(
                        Chem.MolToSmiles(mol, isomericSmiles=True)
                    )
            except Exception as e:
                raise RuntimeError(f"Failed to process product: {e}")

        if not non_cofactor_reactants or not non_cofactor_products:
            return ">>"

        filtered_reaction_smiles = (
            ".".join(non_cofactor_reactants) + ">>" + ".".join(non_cofactor_products)
        )

        new_rxn = AllChem.ReactionFromSmarts(filtered_reaction_smiles, useSmiles=True)
        if new_rxn is None:
            raise ValueError(f"Invalid filtered reaction string: {filtered_reaction_smiles}")

        _clean_up_atom_maps(new_rxn)

        try:
            # Compose a SMIRKS-like output using SMARTS templates so atom maps survive.
            reactant_smarts = [
                Chem.MolToSmarts(m, isomericSmiles=True) for m in new_rxn.GetReactants()
            ]
            product_smarts = [
                Chem.MolToSmarts(m, isomericSmiles=True) for m in new_rxn.GetProducts()
            ]
            return ".".join(reactant_smarts) + ">>" + ".".join(product_smarts)
        except Exception as e:
            raise RuntimeError(f"Failed to convert filtered reaction to SMARTS: {e}")

    except Exception as e:
        raise RuntimeError(f"Failed to remove cofactors from reaction: {rxn_smiles}, Error: {e}")


def contains_cofactor(rxn_smiles: str) -> bool:
    """
    Return True if any reactant or product hashes to a ubiquitous metabolite.
    """
    try:
        native_metabolites = _load_native_metabolites()

        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            raise ValueError(f"Invalid reaction string: {rxn_smiles}")

        for mol in rxn.GetReactants():
            try:
                if get_molecule_hash_from_mol(mol) in native_metabolites:
                    return True
            except Exception as e:
                raise RuntimeError(f"Failed to process reactant: {e}")

        for mol in rxn.GetProducts():
            try:
                if get_molecule_hash_from_mol(mol) in native_metabolites:
                    return True
            except Exception as e:
                raise RuntimeError(f"Failed to process product: {e}")

        return False

    except Exception as e:
        raise RuntimeError(f"Failed to check cofactors in reaction: {rxn_smiles}, Error: {e}")