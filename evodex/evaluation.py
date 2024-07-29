import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from evodex.synthesis import project_evodex_operator
from evodex.formula import calculate_formula_diff
from evodex.utils import get_molecule_hash
import json
from itertools import combinations
import re

# Initialize caches
evodex_f_cache = {}
evodex_data_cache = {}

def _add_hydrogens(smirks):
    """
    Add hydrogens to both sides of the SMIRKS.

    Parameters:
    smirks (str): The SMIRKS string representing the reaction.

    Returns:
    str: The SMIRKS string with hydrogens added to both the substrate and product.
    """
    try:
        substrate, product = smirks.split('>>')
        substrate_mol = Chem.MolFromSmiles(substrate)
        product_mol = Chem.MolFromSmiles(product)
        if substrate_mol is None or product_mol is None:
            raise ValueError(f"Invalid SMILES in SMIRKS: {smirks}")
        substrate_mol = Chem.AddHs(substrate_mol)
        product_mol = Chem.AddHs(product_mol)
        substrate_smiles = Chem.MolToSmiles(substrate_mol)
        product_smiles = Chem.MolToSmiles(product_mol)
        smirks_with_h = f"{substrate_smiles}>>{product_smiles}"
        return smirks_with_h
    except Exception as e:
        print(f"Error in _add_hydrogens: {e}")
        return smirks

def assign_evodex_F(smirks):
    """
    Assign an EVODEX-F ID to a given SMIRKS.

    The function first adds hydrogens to both the substrate and product sides of the SMIRKS. 
    It then calculates the difference in molecular formulas between the substrate and the product.
    This formula difference is used to search a pre-loaded cache of EVODEX-F IDs to find a match.

    Parameters:
    smirks (str): The SMIRKS string representing the reaction.

    Returns:
    str: The EVODEX-F ID, if matched. Returns None if no match is found.

    Example:
    >>> smirks = "CCO>>CC=O"
    >>> assign_evodex_F(smirks)
    """
    smirks_with_h = _add_hydrogens(smirks)
    formula_diff = calculate_formula_diff(smirks_with_h)
    print(f"Formula difference: {formula_diff}")
    evodex_f = _load_evodex_f()
    frozen_f_query = frozenset(formula_diff.items())
    evodex_f_id = evodex_f.get(frozen_f_query)
    print(f"Matched EVODEX-F ID: {evodex_f_id}")
    return evodex_f_id

def _load_evodex_f():
    """
    Load EVODEX-F cache from the CSV file.

    Returns:
    dict: A dictionary mapping formula differences to EVODEX-F IDs.

    Raises:
    FileNotFoundError: If the EVODEX-F CSV file is not found.
    """
    global evodex_f_cache
    if not evodex_f_cache:
        script_dir = os.path.dirname(__file__)
        rel_path = os.path.join('..', 'evodex/data', 'EVODEX-F_unique_formulas.csv')
        filepath = os.path.abspath(os.path.join(script_dir, rel_path))
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        evodex_f_df = pd.read_csv(filepath)
        for index, row in evodex_f_df.iterrows():
            formula_diff = eval(row['formula'])
            evodex_id = row['id']
            sources = _parse_sources(row['sources'])
            if frozenset(formula_diff.items()) not in evodex_f_cache:
                evodex_f_cache[frozenset(formula_diff.items())] = []
            evodex_f_cache[frozenset(formula_diff.items())].append(evodex_id)
    # print(f"Loaded EVODEX-F data: {evodex_f_cache}")
    return evodex_f_cache

def _parse_sources(sources):
    """
    Parse the sources field from the CSV file.

    Parameters:
    sources (str): The sources field as a string.

    Returns:
    list: A list of source strings.
    """
    sources = sources.replace('"', '')  # Remove all double quotes
    return sources.split(',')  # Split by commas

def match_operators(smirks, evodex_type='E'):
    """
    Assign complete-style operators based on a given SMIRKS and EVODEX type.

    This function splits the reaction SMIRKS into substrates and products,
    enumerates all possible pairings, and runs a helper method to find
    matching operators for each pairing.

    Parameters:
    smirks (str): The SMIRKS string representing the reaction.
    evodex_type (str): The type of EVODEX operator (Electronic 'E', Nearest-Neighbor 'N', 
    or Core 'C', default is 'E').

    Returns:
    list: A list of valid operator IDs. Returns an empty list if no matching operators are found.
    """
    # Initialize valid operators list
    valid_operators = []

    try:
        # Split the SMIRKS string into substrates and products
        if '>>' in smirks:
            substrates, products = smirks.split('>>')
            substrate_list = substrates.split('.')
            product_list = products.split('.')

            # Assign an integer index to each substrate and product
            substrate_indices = list(range(len(substrate_list)))
            product_indices = list(range(len(product_list)))

            # Construct new reaction objects combinatorially
            all_pairings = set()
            for i in range(1, len(substrate_indices) + 1):
                for j in range(1, len(product_indices) + 1):
                    for reactant_combo in combinations(substrate_indices, i):
                        for product_combo in combinations(product_indices, j):
                            all_pairings.add((frozenset(reactant_combo), frozenset(product_combo)))

            # Generate all pairings of substrates and products
            for pairing in all_pairings:
                reactant_indices, product_indices = pairing
                reactant_smiles = '.'.join([substrate_list[i] for i in sorted(reactant_indices)])
                product_smiles = '.'.join([product_list[i] for i in sorted(product_indices)])
                pairing_smirks = f"{reactant_smiles}>>{product_smiles}"
                valid_operators.extend(_match_operator(pairing_smirks, evodex_type))

    except Exception as e:
        print(f"Error processing SMIRKS {smirks}: {e}")
    
    return valid_operators

def _match_operator(smirks, evodex_type='E'):
    """
    Helper function to assign a complete-style operator based on a given SMIRKS and EVODEX type.

    Parameters:
    smirks (str): The SMIRKS string representing the reaction.
    evodex_type (str): The type of EVODEX operator (Electronic 'E', Nearest-Neighbor 'N', 
    or Core 'C', default is 'E').

    Returns:
    list: A list of valid operator IDs. Returns an empty list if no matching operators are found.
    """
    # Calculate the formula difference
    smirks_with_h = _add_hydrogens(smirks)
    formula_diff = calculate_formula_diff(smirks_with_h)
    print(f"Formula difference for {smirks}: {formula_diff}")

    # Lazy load the operators associated with each formula
    evodex_f = _load_evodex_f()
    if evodex_f is None:
        return []

    f_id_list = evodex_f.get(frozenset(formula_diff.items()), [])
    if not f_id_list:
        print(f"No F IDs found for formula {formula_diff}")
        return []
    f_id = f_id_list[0]  # Extract the first F_id from the list

    evodex_data = _load_evodex_data()

    if f_id not in evodex_data:
        print(f"F ID {f_id} not in EVODEX data")
        return []

    # Retrieve all operators of the right type associated with the formula difference
    potential_operators = evodex_data[f_id].get(evodex_type, [])
    evodex_ids = [op["id"] for op in potential_operators]
    print(f"Potential operator IDs for {smirks} of type {evodex_type}: {evodex_ids}")

    # Split the input smirks into substrates and products
    sub_smiles, pdt_smiles = smirks.split('>>')

    # Convert pdt_smiles to a hash
    pdt_hash = get_molecule_hash(pdt_smiles)

    # Iterate through potential operators and test
    valid_operators = []
    for operator in potential_operators:
        try:
            id = operator["id"]
            projected_pdts = project_evodex_operator(id, sub_smiles)
            for proj_smiles in projected_pdts:
                proj_hash = get_molecule_hash(proj_smiles)
                if proj_hash == pdt_hash:
                    valid_operators.append(id)
        except Exception as e:
            print(f"{operator['id']} errored: {e}")

    return valid_operators

def _load_evodex_data():
    """
    Return pre-cached JSON object.

    Returns:
    dict: The loaded EVODEX data as a dictionary.

    Raises:
    FileNotFoundError: If the JSON file is not found.
    """
    global evodex_data_cache
    if evodex_data_cache is not None:
        print("Using cached EVODEX data.")
        return evodex_data_cache

    # Load the EVODEX data from the JSON file and return it as an object.
    script_dir = os.path.dirname(__file__)
    rel_path = os.path.join('..', 'evodex/data', 'evaluation_operator_data.json')
    json_filepath = os.path.abspath(os.path.join(script_dir, rel_path))
    if os.path.exists(json_filepath):
        print(f"Loading EVODEX data from {json_filepath}.")
        with open(json_filepath, 'r') as json_file:
            evodex_data_cache = json.load(json_file)
        return evodex_data_cache
    
    # If the JSON file does not exist, create the data from CSV files and save it as JSON.
    print(f"{json_filepath} does not exist. Creating new JSON data from CSV files.")
    
    # Index EVODEX data as JSON files
    e_data = _create_evodex_json('E')
    n_data = _create_evodex_json('N')
    c_data = _create_evodex_json('C')

    # Initialize cache
    evodex_data_cache = {}

    # Load EVODEX-F data
    rel_path = os.path.join('..', 'evodex/data', 'EVODEX-F_unique_formulas.csv')
    csv_filepath = os.path.abspath(os.path.join(script_dir, rel_path))
    evodex_f_df = pd.read_csv(csv_filepath)

    for index, row in evodex_f_df.iterrows():
        f_id = row['id']
        p_ids = _parse_sources(row['sources'])
        all_operator_data_for_F_line = {"C": [], "N": [], "E": []}

        for p_id in p_ids:
            if p_id in c_data and c_data[p_id] not in all_operator_data_for_F_line["C"]:
                all_operator_data_for_F_line["C"].append(c_data[p_id])
            if p_id in n_data and n_data[p_id] not in all_operator_data_for_F_line["N"]:
                all_operator_data_for_F_line["N"].append(n_data[p_id])
            if p_id in e_data and e_data[p_id] not in all_operator_data_for_F_line["E"]:
                all_operator_data_for_F_line["E"].append(e_data[p_id])

        evodex_data_cache[f_id] = all_operator_data_for_F_line

    # Save the combined EVODEX data to a JSON file
    with open(json_filepath, 'w') as json_file:
        json.dump(evodex_data_cache, json_file, indent=4)
    print(f"Saved EVODEX data to {json_filepath}.")

    return evodex_data_cache

def _create_evodex_json(file_suffix):
    """
    Create a dictionary from EVODEX CSV files and save as JSON.

    Parameters:
    file_suffix (str): The suffix for the EVODEX data type ('E', 'N', or 'C').

    Returns:
    dict: The created dictionary from the EVODEX CSV files.

    Raises:
    FileNotFoundError: If the CSV file is not found.
    """
    script_dir = os.path.dirname(__file__)
    csv_path = os.path.join('..', f'evodex/data/EVODEX-{file_suffix}_reaction_operators.csv')
    json_path = os.path.join('..', f'evodex/data/evodex_{file_suffix.lower()}_data.json')

    csv_filepath = os.path.abspath(os.path.join(script_dir, csv_path))
    json_filepath = os.path.abspath(os.path.join(script_dir, json_path))

    if not os.path.exists(csv_filepath):
        raise FileNotFoundError(f"File not found: {csv_filepath}")

    evodex_df = pd.read_csv(csv_filepath)

    evodex_dict = {}
    for index, row in evodex_df.iterrows():
        evodex_id = row['id']
        sources = _parse_sources(row['sources'])
        for source in sources:
            if source not in evodex_dict:
                evodex_dict[source] = []
            evodex_dict[source].append({
                "id": evodex_id,
                "smirks": row['smirks']
            })

    with open(json_filepath, 'w') as json_file:
        json.dump(evodex_dict, json_file, indent=4)

    print(f"EVODEX-{file_suffix} data has been saved to {json_filepath}")
    return evodex_dict

# Example usage:
if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    smirks = "CCCO>>CCC=O"
    is_valid_formula = assign_evodex_F(smirks)
    print(f"{smirks} matches: {is_valid_formula}")

    matching_operators = match_operators(smirks, 'C')
    print(f"Matching operators for {smirks}: {matching_operators}")
