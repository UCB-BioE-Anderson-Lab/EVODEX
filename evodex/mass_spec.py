import os
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pandas as pd
import json
from evodex.synthesis import project_evodex_operator
from evodex.evaluation import _load_evodex_data, _parse_sources
from evodex.utils import get_molecule_hash

# Initialize caches
evodex_m_cache = None
evodex_data_cache = None
evodex_m_to_f_cache = None

def calculate_mass(smiles):
    """
    Calculate the molecular mass of a given SMILES string.

    Parameters:
    smiles (str): The SMILES string representing the molecule.

    Returns:
    float: The exact mass of the molecule.

    Raises:
    ValueError: If the SMILES string is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return rdMolDescriptors.CalcExactMolWt(mol)
    else:
        raise ValueError(f"Invalid SMILES string: {smiles}")

def _load_evodex_m():
    """
    Load EVODEX-M cache from the CSV file.

    Returns:
    list: A list of dictionaries, each containing 'id', 'mass', and 'sources' keys.

    Raises:
    FileNotFoundError: If the EVODEX-M CSV file is not found.
    """
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
    """
    Load or create EVODEX-M to EVODEX-F mapping.

    Returns:
    dict: A dictionary mapping EVODEX-M IDs to lists of EVODEX-F IDs.

    Raises:
    FileNotFoundError: If the necessary CSV files are not found.
    """
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
    """
    Find EVODEX-M entries that correspond to a given mass difference within a specified precision.

    Parameters:
    mass_diff (float): The mass difference to search for.
    precision (float): The precision within which to match the mass difference (default is 0.01).

    Returns:
    list: A list of matching EVODEX-M entries, each containing 'id' and 'mass'.
    """
    evodex_m = _load_evodex_m()
    matching_entries = [
        {"id": entry["id"], "mass": entry["mass"]}
        for entry in evodex_m
        if abs(entry["mass"] - mass_diff) <= precision
    ]
    return matching_entries

def get_reaction_operators(mass_diff, precision=0.01):
    """
    Retrieve reaction operators that could explain the mass difference.

    Parameters:
    mass_diff (float): The mass difference to search for.
    precision (float): The precision within which to match the mass difference (default is 0.01).

    Returns:
    tuple: A dictionary of matching reaction operators by type ('E', 'C', 'N'), 
           a list of matching EVODEX-M entries, and a list of corresponding EVODEX-F IDs.
    """
    matching_evodex_m = find_evodex_m(mass_diff, precision)
    if not matching_evodex_m:
        return {}, matching_evodex_m, []

    evodex_m_to_f = _load_evodex_m_to_f()
    evodex_data = _load_evodex_data()

    matching_operators = {"E": [], "C": [], "N": []}
    f_ids = set()
    for entry in matching_evodex_m:
        evodex_m_id = entry["id"]
        if evodex_m_id in evodex_m_to_f:
            for f_id in evodex_m_to_f[evodex_m_id]:
                f_ids.add(f_id)
                if f_id in evodex_data:
                    for op_type, ops in evodex_data[f_id].items():
                        matching_operators[op_type].extend(ops)

    return matching_operators, matching_evodex_m, list(f_ids)


def predict_products(smiles, mass_diff, precision=0.01, evodex_type='E'):
    """
    Project all EVODEX-E, -N, or -C operators consistent with the EVODEX-M onto given substrates and predict the products.

    Parameters:
    smiles (str): The SMILES string representing the substrate.
    mass_diff (float): The mass difference to search for.
    precision (float): The precision within which to match the mass difference (default is 0.01).
    evodex_type (str): The type of EVODEX operator (Electronic 'E', Nearest-Neighbor 'N', 
    or Core 'C', default is 'E').

    Returns:
    dict: A dictionary of predicted products with details of the projections.

    Example:
    >>> smiles = "CCCO"
    >>> mass_diff = 14.016
    >>> precision = 0.01
    >>> predict_products(smiles, mass_diff, precision, 'E')
    """    
    matching_operators, matching_evodex_m, f_ids = get_reaction_operators(mass_diff, precision)
    evodex_e_ops = matching_operators.get(evodex_type, [])

    results = {}
    for operator in evodex_e_ops:
        try:
            # Accessing the first element in the matching_evodex_m list
            m_id = matching_evodex_m[0]["id"]
            operator_id = operator["id"]
            projected_pdts = project_evodex_operator(operator_id, smiles)
            for proj_smiles in projected_pdts:
                proj_hash = get_molecule_hash(proj_smiles)
                if proj_hash not in results:
                    results[proj_hash] = {
                        "smiles": proj_smiles,
                        "projections": {}
                    }
                formula_mass_key = (tuple(f_ids), m_id)
                if formula_mass_key not in results[proj_hash]["projections"]:
                    results[proj_hash]["projections"][formula_mass_key] = []
                if operator_id not in results[proj_hash]["projections"][formula_mass_key]:
                    results[proj_hash]["projections"][formula_mass_key].append(operator_id)
        except Exception as e:
            print(f"{operator['id']} errored: {str(e)}")

    return results

if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    substrate = "CCCO"
    mass_diff = 14.016  # Example mass difference
    precision = 0.01

    evodex_m = find_evodex_m(mass_diff, precision)
    print(f"Found matching {mass_diff}: {evodex_m}")

    matching_operators, _, _ = get_reaction_operators(mass_diff, precision)
    print(f"Matching operators for mass difference {mass_diff}: {[op['id'] for op_list in matching_operators.values() for op in op_list]}")

    results = predict_products(substrate, mass_diff, precision, 'E')
    for product, details in results.items():
        print(f"Product: {details['smiles']}")
        for (f_ids, m_id), operators in details['projections'].items():
            print(f"  EVODEX-F (formula) IDs: {f_ids}, EVODEX-M (mass) ID: {m_id}, Operators: {operators}")
