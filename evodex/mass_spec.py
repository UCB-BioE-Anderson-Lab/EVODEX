import os
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pandas as pd
import json
from evodex.synthesis import project_evodex_operator
from evodex.evaluation import _load_evodex_data, _parse_sources

# Initialize caches
evodex_m_cache = None
evodex_data_cache = None
evodex_m_to_f_cache = None

def calculate_mass(smiles):
    """Calculate the molecular mass of a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return rdMolDescriptors.CalcExactMolWt(mol)
    else:
        raise ValueError(f"Invalid SMILES string: {smiles}")

def _load_evodex_m():
    """Load EVODEX-M cache from the CSV file."""
    global evodex_m_cache
    if evodex_m_cache is None:
        evodex_m_cache = []
        script_dir = os.path.dirname(__file__)
        rel_path = os.path.join('..', 'evodex/data', 'EVODEX-M_mass_spec_subset.csv')
        filepath = os.path.abspath(os.path.join(script_dir, rel_path))
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        evodex_m_df = pd.read_csv(filepath)
        for index, row in evodex_m_df.iterrows():
            evodex_m_cache.append({
                "id": row['id'],
                "mass": row['mass'],
                "sources": _parse_sources(row['sources'])
            })
    return evodex_m_cache

def _load_evodex_m_to_f():
    """Load or create EVODEX-M to EVODEX-F mapping."""
    global evodex_m_to_f_cache
    if evodex_m_to_f_cache is not None:
        return evodex_m_to_f_cache

    script_dir = os.path.dirname(__file__)
    rel_path = os.path.join('..', 'evodex/data', 'evodex_m_to_F_mapping.csv')  # Correct file path for saving
    filepath = os.path.abspath(os.path.join(script_dir, rel_path))

    if os.path.exists(filepath):
        evodex_m_to_f_cache = {}
        evodex_m_to_f_df = pd.read_csv(filepath)
        for index, row in evodex_m_to_f_df.iterrows():
            if row['evodex_m'] not in evodex_m_to_f_cache:
                evodex_m_to_f_cache[row['evodex_m']] = []
            evodex_m_to_f_cache[row['evodex_m']].append(row['evodex_f'])
        return evodex_m_to_f_cache

    # Load EVODEX-M data
    evodex_m_cache = _load_evodex_m()
    
    # Create the EVODEX-P to EVODEX-F mapping
    script_dir = os.path.dirname(__file__)
    rel_path = os.path.join('..', 'evodex/data', 'EVODEX-F_unique_formulas.csv')
    filepath_f = os.path.abspath(os.path.join(script_dir, rel_path))
    if not os.path.exists(filepath_f):
        raise FileNotFoundError(f"File not found: {filepath_f}")

    evodex_f_df = pd.read_csv(filepath_f)
    p_to_f_map = {}
    for index, row in evodex_f_df.iterrows():
        f_id = row['id']
        p_ids = _parse_sources(row['sources'])
        for p_id in p_ids:
            if p_id not in p_to_f_map:
                p_to_f_map[p_id] = []
            p_to_f_map[p_id].append(f_id)

    # Create the EVODEX-M to EVODEX-F mapping
    evodex_m_to_f_cache = {}
    for entry in evodex_m_cache:
        evodex_m_id = entry["id"]
        for p_id in entry["sources"]:
            if p_id in p_to_f_map:
                if evodex_m_id not in evodex_m_to_f_cache:
                    evodex_m_to_f_cache[evodex_m_id] = []
                evodex_m_to_f_cache[evodex_m_id].extend(p_to_f_map[p_id])

    # Save the mapping to a CSV file
    with open(filepath, 'w') as f:  # Use the correct filepath for saving
        f.write("evodex_m,evodex_f\n")
        for evodex_m_id, evodex_f_ids in evodex_m_to_f_cache.items():
            for evodex_f_id in evodex_f_ids:
                f.write(f"{evodex_m_id},{evodex_f_id}\n")

    return evodex_m_to_f_cache

def find_evodex_m(mass_diff, precision=0.01):
    """Find EVODEX-M entries that correspond to a given mass difference within a specified precision."""
    evodex_m = _load_evodex_m()
    matching_entries = [
        {"id": entry["id"], "mass": entry["mass"]}
        for entry in evodex_m
        if abs(entry["mass"] - mass_diff) <= precision
    ]
    return matching_entries

def get_reaction_operators(mass_diff, precision=0.01):
    """Retrieve reaction operators that could explain the mass difference."""
    matching_evodex_m = find_evodex_m(mass_diff, precision)
    if not matching_evodex_m:
        return {}

    evodex_m_to_f = _load_evodex_m_to_f()
    evodex_data = _load_evodex_data()

    matching_operators = {"E": [], "C": [], "N": []}
    for entry in matching_evodex_m:
        evodex_m_id = entry["id"]
        if evodex_m_id in evodex_m_to_f:
            f_ids = evodex_m_to_f[evodex_m_id]
            for f_id in f_ids:
                if f_id in evodex_data:
                    for op_type, ops in evodex_data[f_id].items():
                        matching_operators[op_type].extend(ops)

    return matching_operators

def project_mass_diff_operators(smiles, mass_diff, precision=0.01):
    """Project all EVODEX-E operators consistent with the EVODEX-M onto given substrates and predict the products."""
    matching_operators = get_reaction_operators(mass_diff, precision)
    evodex_e_ops = matching_operators.get("E", [])

    valid_products = []
    for operator in evodex_e_ops:
        try:
            id = operator["id"]
            projected_pdts = project_evodex_operator(id, smiles)
            valid_products.extend(projected_pdts)
        except Exception as e:
            print(f"{operator['id']} errored: {str(e)}")

    return valid_products

# Example usage
if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    substrate = "CCCO"
    mass_diff = 14.016  # Example mass difference
    precision = 0.01

    evodex_m = find_evodex_m(mass_diff, precision)
    print(f"Found matching {mass_diff}: {evodex_m}")

    matching_operators = get_reaction_operators(mass_diff, precision)
    print(f"Matching operators for mass difference {mass_diff}: {[op['id'] for op_list in matching_operators.values() for op in op_list]}")

    predicted_products = project_mass_diff_operators(substrate, mass_diff, precision)
    print(f"Predicted products for substrate {substrate} and mass difference {mass_diff}: {predicted_products}")
