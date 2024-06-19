import csv
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Set

# Global variable to store native metabolites set
_native_metabolites = None

def _load_native_metabolites() -> Set[str]:
    global _native_metabolites
    if _native_metabolites is not None:
        return _native_metabolites

    native_metabolites_file = 'data/2024_06_18-Native_Metabolites.tsv'
    _native_metabolites = set()
    try:
        with open(native_metabolites_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                inchi = row['inchi'].strip().strip('"')
                mol = Chem.MolFromInchi(inchi)
                if mol:
                    canonical_smiles = canonicalize_molecule(mol)
                    _native_metabolites.add(canonical_smiles)
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

def canonicalize_molecule(mol: Chem.Mol) -> str:
    try:
        mol_copy = Chem.Mol(mol)  # Make a copy of the molecule
        for atom in mol_copy.GetAtoms():  # Clear atom maps
            atom.SetAtomMapNum(0)
        for atom in mol_copy.GetAtoms():  # Update implicit valence
            atom.UpdatePropertyCache()
        mol_copy = Chem.AddHs(mol_copy)  # Add explicit hydrogens
        Chem.SanitizeMol(mol_copy)  # Sanitize the molecule
        canonical_smiles = Chem.MolToSmiles(mol_copy, canonical=True)
    except Exception as e:
        raise RuntimeError(f"Failed to canonicalize molecule: {e}")
    return canonical_smiles

def remove_cofactors(smirks: str) -> str:
    try:
        native_metabolites = _load_native_metabolites()

        # Load the input SMIRKS as a reaction object
        rxn = AllChem.ReactionFromSmarts(smirks)
        if not rxn:
            raise ValueError(f"Invalid SMIRKS string: {smirks}")

        # Identify non-cofactor reactants and products
        non_cofactor_reactants = []
        non_cofactor_products = []

        for mol in rxn.GetReactants():
            try:
                canonical_smiles = canonicalize_molecule(mol)
                if canonical_smiles and canonical_smiles not in native_metabolites:
                    non_cofactor_reactants.append(Chem.MolToSmiles(mol, isomericSmiles=True))
            except Exception as e:
                raise RuntimeError(f"Failed to process reactant: {e}")

        for mol in rxn.GetProducts():
            try:
                canonical_smiles = canonicalize_molecule(mol)
                if canonical_smiles and canonical_smiles not in native_metabolites:
                    non_cofactor_products.append(Chem.MolToSmiles(mol, isomericSmiles=True))
            except Exception as e:
                raise RuntimeError(f"Failed to process product: {e}")

        if not non_cofactor_reactants or not non_cofactor_products:
            raise ValueError("No valid non-cofactor reactants or products found")

        # Create a new reaction with non-cofactor molecules
        reactant_smiles = '.'.join(non_cofactor_reactants)
        product_smiles = '.'.join(non_cofactor_products)
        new_reaction_smirks = f"{reactant_smiles}>>{product_smiles}"

        # Process the new reaction to clean up atom maps
        new_rxn = AllChem.ReactionFromSmarts(new_reaction_smirks)
        if not new_rxn:
            raise ValueError(f"Invalid new reaction SMIRKS: {new_reaction_smirks}")

        _clean_up_atom_maps(new_rxn)
        return AllChem.ReactionToSmarts(new_rxn)
    except Exception as e:
        raise RuntimeError(f"Failed to remove cofactors from SMIRKS: {smirks}, Error: {e}")
