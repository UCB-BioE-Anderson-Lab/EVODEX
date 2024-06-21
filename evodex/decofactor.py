import csv
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Set
from evodex.astatine import hydrogen_to_astatine_molecule
from evodex.utils import calculate_inchi

# Global variable to store native metabolites set
_native_metabolites = None

def _load_native_metabolites() -> Set[str]:
    global _native_metabolites
    if _native_metabolites is not None:
        return _native_metabolites

    script_dir = os.path.dirname(__file__)
    native_metabolites_file = os.path.join(script_dir, 'data', '2024_06_18-Native_Metabolites.tsv')
    _native_metabolites = set()
    try:
        with open(native_metabolites_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                inchi = row['inchi'].strip().strip('"')
                mol = Chem.MolFromInchi(inchi)
                if mol:
                    astatine_mol = hydrogen_to_astatine_molecule(mol)
                    astatine_inchi = Chem.MolToInchiKey(astatine_mol)
                    _native_metabolites.add(astatine_inchi)
    except Exception as e:
        raise RuntimeError(f"Failed to load native metabolites: {e}")
    return _native_metabolites

def _clean_up_atom_maps(rxn: AllChem.ChemicalReaction):
    try:
        substrate_atom_maps = set()

        # Collect atom maps from reactants
        for mol in rxn.GetReactants():
            for atom in mol.GetAtoms():
                atom_map_num = atom.GetAtomMapNum()
                if atom_map_num > 0:
                    substrate_atom_maps.add(atom_map_num)

        # Adjust atom maps in products
        for mol in rxn.GetProducts():
            for atom in mol.GetAtoms():
                atom_map_num = atom.GetAtomMapNum()
                if atom_map_num > 0:
                    if atom_map_num not in substrate_atom_maps:
                        atom.SetAtomMapNum(0)
                    else:
                        substrate_atom_maps.remove(atom_map_num)

        # Adjust atom maps in reactants
        for mol in rxn.GetReactants():
            for atom in mol.GetAtoms():
                atom_map_num = atom.GetAtomMapNum()
                if atom_map_num in substrate_atom_maps:
                    atom.SetAtomMapNum(0)
    except Exception as e:
        raise RuntimeError(f"Failed to clean up atom maps: {e}")

def remove_cofactors(smiles: str) -> str:
    try:
        native_metabolites = _load_native_metabolites()

        # Load the input SMILES as a reaction object
        rxn = AllChem.ReactionFromSmarts(smiles, useSmiles=True)
        if not rxn:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # Identify non-cofactor reactants and products
        non_cofactor_reactants = []
        non_cofactor_products = []

        for mol in rxn.GetReactants():
            try:
                inchi = Chem.MolToInchiKey(mol)
                if inchi and inchi not in native_metabolites:
                    non_cofactor_reactants.append(Chem.MolToSmiles(mol, isomericSmiles=True))
            except Exception as e:
                raise RuntimeError(f"Failed to process reactant: {e}")

        for mol in rxn.GetProducts():
            try:
                inchi = Chem.MolToInchiKey(mol)
                if inchi and inchi not in native_metabolites:
                    non_cofactor_products.append(Chem.MolToSmiles(mol, isomericSmiles=True))
            except Exception as e:
                raise RuntimeError(f"Failed to process product: {e}")

        if not non_cofactor_reactants or not non_cofactor_products:
            return ">>"  # Return an empty smiles when no valid non-cofactor reactants or products are found

        # Create a new reaction with non-cofactor molecules
        reactant_smiles = '.'.join(non_cofactor_reactants)
        product_smiles = '.'.join(non_cofactor_products)
        new_reaction_smiles = f"{reactant_smiles}>>{product_smiles}"

        # Process the new reaction to clean up atom maps
        new_rxn = AllChem.ReactionFromSmarts(new_reaction_smiles, useSmiles=True)
        if not new_rxn:
            raise ValueError(f"Invalid new reaction SMILES: {new_reaction_smiles}")

        _clean_up_atom_maps(new_rxn)

        try:
            reaction_smarts = AllChem.ReactionToSmarts(new_rxn)
            reactant_smiles = [Chem.MolToSmarts(mol, isomericSmiles=True) for mol in new_rxn.GetReactants()]
            product_smiles = [Chem.MolToSmarts(mol, isomericSmiles=True) for mol in new_rxn.GetProducts()]
            modified_smiles = '>>'.join(['.'.join(reactant_smiles), '.'.join(product_smiles)])
        except Exception as e:
            raise ValueError(f"Error converting modified reaction to SMIRKS: {reaction_smarts}") from e

        return modified_smiles
    except Exception as e:
        raise RuntimeError(f"Failed to remove cofactors from SMILES: {smiles}, Error: {e}")
