import csv
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import inchi
from typing import Set
from evodex.utils import get_molecule_hash

# Global variable to store native metabolites set
_native_metabolites = None

def _load_native_metabolites() -> Set[str]:
    global _native_metabolites
    if _native_metabolites is not None:
        # print("Native metabolites already loaded.")
        return _native_metabolites

    # print("----> loading native metabolites from file")
    script_dir = os.path.dirname(__file__)
    ubiquitous_metabolites_file = os.path.join(script_dir, 'data', 'ubiquitous_metabolites.txt')
    backup_file = os.path.join(script_dir, 'data', 'ubiquitous_metabolites_backup.txt')
    _native_metabolites = set()
    
    try:
        with open(ubiquitous_metabolites_file, 'r') as file:
            # print("reading ubiquitous metabolites")
            reader = csv.DictReader(file, delimiter='\t')
            
            with open(backup_file, 'w', newline='') as backup:
                backup_writer = csv.writer(backup, delimiter='\t')
                backup_writer.writerow(['name', 'original_inchi', 'custom_hash'])
                
                for row in reader:
                    name = row['name'].strip().strip('"')
                    # print("processing: ", name)
                    inchi_str = row['inchi'].strip().strip('"')
                    mol = inchi.MolFromInchi(inchi_str)
                    
                    if mol:
                        try:
                            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                            custom_hash = get_molecule_hash(smiles)
                            # print("hash:", custom_hash)
                            _native_metabolites.add(custom_hash)
                            backup_writer.writerow([name, inchi_str, custom_hash])
                        except Exception as e:
                            print(f"Failed to process molecule: {name}, Error: {e}")
    
    except Exception as e:
        raise RuntimeError(f"Failed to load ubiquitous metabolites: {e}")
    
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
    # print("remove cofactors invoked")
    try:
        native_metabolites = _load_native_metabolites()
        # print("post loading metabolites")

        # Load the input SMILES as a reaction object
        rxn = AllChem.ReactionFromSmarts(smiles, useSmiles=True)
        if not rxn:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # Identify non-cofactor reactants and products
        non_cofactor_reactants = []
        non_cofactor_products = []

        # print("----> hashing reactants and products")
        for mol in rxn.GetReactants():
            try:
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                custom_hash = get_molecule_hash(smiles)
                # print("\ntesting reactant hash: ", custom_hash)
                if custom_hash and custom_hash not in native_metabolites:
                    retained_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                    # print("Retained smiles:", retained_smiles)
                    non_cofactor_reactants.append(retained_smiles)       
            except Exception as e:
                raise RuntimeError(f"Failed to process reactant: {e}")

        for mol in rxn.GetProducts():
            try:
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                custom_hash = get_molecule_hash(smiles)
                # print("\ntesting product hash: ", custom_hash)
                if custom_hash and custom_hash not in native_metabolites:
                    retained_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                    # print("Retained smiles:", retained_smiles)
                    non_cofactor_products.append(retained_smiles)
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