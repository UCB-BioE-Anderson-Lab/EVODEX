# evodex/predictions/synthesis.py

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Global cache for EVODEX data
evodex_data_cache = {}

def project_reaction_operator(smirks, substrate):
    """Apply a reaction operator (SMIRKS) to a substrate."""
    # Create a reaction object from the SMIRKS string
    rxn = AllChem.ReactionFromSmarts(smirks)
    if not rxn:
        raise ValueError(f"Failed to create reaction from SMIRKS: {smirks}")
    
    # Check if the reaction has valid reactants and products
    if rxn.GetNumReactantTemplates() == 0 or rxn.GetNumProductTemplates() == 0:
        raise ValueError(f"Reaction has no valid reactants or products: {smirks}")

    # Create a molecule object from the substrate SMILES string
    substrate_mol = Chem.MolFromSmiles(substrate)
    if not substrate_mol:
        raise ValueError(f"Failed to create molecule from substrate SMILES: {substrate}")
    
    # Add hydrogens to the substrate molecule
    substrate_mol = Chem.AddHs(substrate_mol)

    # Apply the reaction to the substrate molecule
    products = rxn.RunReactants((substrate_mol,))
    if not products:
        return []
    
    # Convert the product molecules to canonical SMILES strings
    unique_products = set()
    for product_tuple in products:
        for product in product_tuple:
            if product:
                canonical_smiles = Chem.MolToSmiles(product)
                unique_products.add(canonical_smiles)
    
    return list(unique_products)

def project_evodex_operator(evodex_id, substrate):
    """Apply an EVODEX operator to a substrate."""
    smirks = _lookup_smirks_by_evodex_id(evodex_id)
    if smirks is None:
        raise ValueError(f"SMIRKS not found for EVODEX ID: {evodex_id}")
    return project_reaction_operator(smirks, substrate)

def project_synthesis_operators(substrate):
    """Project all synthesis subset EVODEX-E operators on a substrate."""
    evodex_e_operators = _load_synthesis_evodex_e_operators()
    applicable_operators = {}
    for evodex_id, smirks in evodex_e_operators.items():
        try:
            products = project_reaction_operator(smirks, substrate)
            if products:
                applicable_operators[evodex_id] = products
        except Exception:
            pass
    return applicable_operators

def _lookup_smirks_by_evodex_id(evodex_id):
    """Look up SMIRKS by EVODEX ID."""
    if evodex_id.startswith("E"):
        evodex_df = _load_evodex_data('EVODEX-E_reaction_operators.csv')
    elif evodex_id.startswith("C"):
        evodex_df = _load_evodex_data('EVODEX-C_reaction_operators.csv')
    elif evodex_id.startswith("N"):
        evodex_df = _load_evodex_data('EVODEX-N_reaction_operators.csv')
    else:
        return None

    smirks = evodex_df.loc[evodex_df['id'] == evodex_id, 'smirks'].values
    return smirks[0] if len(smirks) > 0 else None

def _load_synthesis_evodex_e_operators():
    """Load synthesis subset EVODEX-E operators."""
    return _load_evodex_data('EVODEX-E_synthesis_subset.csv', key_column='id', value_column='smirks')

def _load_evodex_data(filename, key_column=None, value_column=None):
    """Lazy-load EVODEX data from a CSV file and cache it."""
    global evodex_data_cache
    if filename not in evodex_data_cache:
        script_dir = os.path.dirname(__file__)
        rel_path = os.path.join('..', 'evodex/data', filename)
        filepath = os.path.abspath(os.path.join(script_dir, rel_path))
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        evodex_df = pd.read_csv(filepath)
        if key_column and value_column:
            evodex_data_cache[filename] = dict(zip(evodex_df[key_column], evodex_df[value_column]))
        else:
            evodex_data_cache[filename] = evodex_df
    return evodex_data_cache[filename]

# Example usage:
if __name__ == "__main__":
    # Run direct projection
    smirks = "[H][C:8]([C:7])([O:9][H])[H:19]>>[C:7][C:8](=[O:9])[H:19]"
    substrate = "CCCO"
    result = project_reaction_operator(smirks, substrate)
    print("Direct projection: ",result)

    # Project from EVODEX ID
    evodex_id = "EVODEX-E170"
    result = project_evodex_operator(evodex_id, substrate)
    print("Refernced EVODEX projection: ",result)

    # Project All Synthesis Subset EVODEX
    result = project_synthesis_operators(substrate)
    print("All Synthesis Subset projection: ",result)
