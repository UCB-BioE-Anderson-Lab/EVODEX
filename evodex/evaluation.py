import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import AllChem
from evodex.synthesis import project_evodex_operator
from evodex.formula import calculate_formula_diff
from evodex.utils import get_molecule_hash
import json
from itertools import combinations
from evodex.projection import project_operator

# Initialize caches
evodex_f_cache = None
evodex_data_cache = None

def _add_hydrogens(smirks):
    """
    Add hydrogens to both sides of the SMIRKS.

    Parameters:
    smirks (str): The SMIRKS string representing the reaction.

    Returns:
    str: The SMIRKS string with hydrogens added to both the substrate and product.
    """
    substrate, product = smirks.split('>>')
    substrate_mol = Chem.MolFromSmiles(substrate)
    product_mol = Chem.MolFromSmiles(product)
    substrate_mol = Chem.AddHs(substrate_mol)
    product_mol = Chem.AddHs(product_mol)
    substrate_smiles = Chem.MolToSmiles(substrate_mol)
    product_smiles = Chem.MolToSmiles(product_mol)
    smirks_with_h = f"{substrate_smiles}>>{product_smiles}"
    return smirks_with_h

def operator_matches_reaction(operator_smirks: str, reaction_smiles: str) -> bool:
    """
    Determine if a reaction operator (SMIRKS) produces the product of a given reaction (SMILES).

    Parameters:
    operator_smirks (str): The reaction operator in SMIRKS format.
    reaction_smiles (str): A full reaction in SMILES format (substrates>>products).

    Returns:
    bool: True if the operator can generate the correct product(s) from the substrate(s), else False.
    """
    # New input validation and logging
    if not isinstance(operator_smirks, str) or not operator_smirks.strip():
        raise ValueError(f"[operator_matches_reaction] Invalid operator_smirks: {operator_smirks}")
    if not isinstance(reaction_smiles, str) or not reaction_smiles.strip():
        raise ValueError(f"[operator_matches_reaction] Invalid reaction_smiles: {reaction_smiles}")

    # String-based parsing and count checks before loading RDKit Mol objects
    op_sub_parts = operator_smirks.split(">>")[0].split(".")
    rxn_sub_parts = reaction_smiles.split(">>")[0].split(".")
    if len(rxn_sub_parts) != len(op_sub_parts):
        return False

    op_prod_parts = operator_smirks.split(">>")[1].split(".")
    rxn_prod_parts = reaction_smiles.split(">>")[1].split(".")
    if len(rxn_prod_parts) != len(op_prod_parts):
        return False

    # Now load operator substrate SMARTS patterns and reaction substrate molecules with explicit Hs
    from itertools import permutations
    op_sub_smarts = [Chem.MolFromSmarts(m) for m in op_sub_parts if Chem.MolFromSmarts(m)]
    rxn_sub_mols = [Chem.AddHs(Chem.MolFromSmiles(m)) for m in rxn_sub_parts if Chem.MolFromSmiles(m)]

    if len(rxn_sub_mols) > len(op_sub_smarts):
        # print("reaction has more substrates than operator expects")
        return False

    op_prod_mols = operator_smirks.split(">>")[1].split(".")
    rxn_prod_mols = reaction_smiles.split(">>")[1].split(".")
    if len(rxn_prod_mols) > len(op_prod_mols):
        # print("reaction has more products than operator expects")
        return False

    # Substructure matching block for substrate parts
    found_match = False
    for op_perm in permutations(op_sub_smarts):
        used_indices = set()
        all_match = True
        for op_mol in op_perm:
            match_found = False
            for i, rxn_mol in enumerate(rxn_sub_mols):
                if i in used_indices:
                    continue
                if rxn_mol.HasSubstructMatch(op_mol):
                    used_indices.add(i)
                    match_found = True
                    break
            if not match_found:
                all_match = False
                break
        if all_match:
            found_match = True
            break

    if not found_match:
        return False

    try:
        reaction_smiles_with_h = _add_hydrogens(reaction_smiles)
        # Parse the reaction, remove atom maps, and convert back to SMILES
        rxn = AllChem.ReactionFromSmarts(reaction_smiles_with_h, useSmiles=True)
        for i in range(rxn.GetNumReactantTemplates()):
            for atom in rxn.GetReactantTemplate(i).GetAtoms():
                atom.SetAtomMapNum(0)
        for i in range(rxn.GetNumProductTemplates()):
            for atom in rxn.GetProductTemplate(i).GetAtoms():
                atom.SetAtomMapNum(0)
        substrates = '.'.join([Chem.MolToSmiles(rxn.GetReactantTemplate(i)) for i in range(rxn.GetNumReactantTemplates())])
        expected_products = '.'.join([Chem.MolToSmiles(rxn.GetProductTemplate(i)) for i in range(rxn.GetNumProductTemplates())])
        expected_hashes = set(get_molecule_hash(p) for p in expected_products.split('.'))

        projected_product_sets = project_operator(operator_smirks, substrates)
        for product_string in projected_product_sets:
            projected_hashes = set(get_molecule_hash(p) for p in product_string.split('.'))
            if projected_hashes == expected_hashes:
                return True
        # print("projections don't match")
        return False
    except Exception as e:
        # print("error")
        # print(e)
        return False
    
def assign_evodex_F(smiles):
    """
    Assign an EVODEX-F ID to a given reaction SMILES.

    The function first adds hydrogens to both the substrate and product sides of the SMILES. 
    It then calculates the difference in molecular formulas between the substrate and the product.
    This formula difference is used to search a pre-loaded cache of EVODEX-F IDs to find a match.

    Parameters:
    smiles (str): The SMILES string representing the reaction.

    Returns:
    str: The EVODEX-F ID, if matched. Returns None if no match is found.

    Example:
    >>> smiles = "CCO>>CC=O"
    >>> assign_evodex_F(smiles)
    """
    smirks_with_h = _add_hydrogens(smiles)
    formula_diff = calculate_formula_diff(smirks_with_h)
    # print("Formula difference:", formula_diff)
    evodex_f = _load_evodex_f()
    evodex_f_id = evodex_f.get(frozenset(formula_diff.items()))
    # print("Matched EVODEX-F ID:", evodex_f_id)
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
    if evodex_f_cache is None:
        evodex_f_cache = {}
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

def match_operators(smiles, evodex_type='E'):
    """
    Assign complete-style operators based on a given reaction SMILES and EVODEX type.

    This function splits the reaction SMILES into substrates and products,
    enumerates all possible pairings, and runs a helper method to find
    matching operators for each pairing.

    Parameters:
    smiles (str): The SMILES string representing the reaction, e.g. 'CCO>>CC=O'
    evodex_type (str): The type of EVODEX operator (Electronic 'E', Nearest-Neighbor 'N', 
    or Core 'C', default is 'E').

    Returns:
    list: A list of valid operator IDs. Returns an empty list if no matching operators are found.
    """
    # Initialize valid operators list
    valid_operators = []

    try:
        # Split the SMILES string into substrates and products
        if '>>' in smiles:
            substrates, products = smiles.split('>>')
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
                pairing_smiles = f"{reactant_smiles}>>{product_smiles}"
                valid_operators.extend(_match_operator(pairing_smiles, evodex_type))

    except Exception as e:
        print(f"Error processing SMILES {smiles}: {e}")
    
    return valid_operators

def _match_operator(smiles, evodex_type='E'):
    """
    Helper function to assign a complete-style operator based on a given SMILES and EVODEX type.

    Parameters:
    smiles (str): The SMILES string representing the reaction.
    evodex_type (str): The type of EVODEX operator (Electronic 'E', Nearest-Neighbor 'N', 
    or Core 'C', default is 'E').

    Returns:
    list: A list of valid operator IDs. Returns an empty list if no matching operators are found.
    """
    # Calculate the formula difference
    smiles_with_h = _add_hydrogens(smiles)
    formula_diff = calculate_formula_diff(smiles_with_h)
    # print("Formula difference:", formula_diff)

    # Lazy load the operators associated with each formula
    evodex_f = _load_evodex_f()
    if evodex_f is None:
        return {}

    f_id_list = evodex_f.get(frozenset(formula_diff.items()), [])
    if not f_id_list:
        return {}
    f_id = f_id_list[0]  # Extract the single F_id from the list

    # print(f"Potential F ID for formula {formula_diff}: {f_id}")

    evodex_data = _load_evodex_data()

    if f_id not in evodex_data:
        return {}

    # Retrieve all operators of the right type associated with the formula difference
    potential_operators = evodex_data[f_id].get(evodex_type, [])
    evodex_ids = [op["id"] for op in potential_operators]
    # print(f"Potential operator IDs for {smiles} of type {evodex_type}: {evodex_ids}")

    # Split the input smiles into substrates and products
    sub_smiles, pdt_smiles = smiles.split('>>')

    # Convert pdt_smiles to a hash
    pdt_hash = get_molecule_hash(pdt_smiles)
    # print(f"Expected product: {pdt_smiles} with hash {pdt_hash}")

    # Iterate through potential operators and test
    valid_operators = []
    for operator in potential_operators:
        try:
            id = operator["id"]
            # print(f"Projecting: {id} on {sub_smiles}")
            projected_pdts = project_operator(operator['smirks'], sub_smiles)
            for proj_smiles in projected_pdts:
                proj_hash = get_molecule_hash(proj_smiles)
                # print(f"Projected product: {proj_smiles} with hash {proj_hash}")
                if proj_hash == pdt_hash:
                    # print("MATCH FOUND!")
                    valid_operators.append(id)
        except Exception as e:
            # print(f"{operator['id']} errored: {e}")
            pass

    return valid_operators


# New function: find_exact_matching_operators
def find_exact_matching_operators(p_smiles, evodex_type='E'):
    """
    Find operators that exactly match the given P reaction.

    Uses formula diff lookup to limit candidates.

    Parameters:
    p_smiles (str): SMILES string of P reaction (substrates>>products)
    evodex_type (str): Operator type ('E', 'N', 'C')

    Returns:
    list: List of operator IDs that exactly match P
    """
    f_id_list = assign_evodex_F(p_smiles)
    if not f_id_list:
        return []
    f_id = f_id_list[0]
    # print(f"Matched F ID: {f_id} for P SMILES: {p_smiles}")

    evodex_data = _load_evodex_data()
    if f_id not in evodex_data:
        return []

    candidate_operators = evodex_data[f_id].get(evodex_type, [])
    substrates, products = p_smiles.split('>>')
    products_hash = get_molecule_hash(products)

    exact_matches = []
    for operator in candidate_operators:
        op_id = operator['id']
        op_smirks = operator.get('smirks', 'MISSING')
        # print(f"Testing operator {op_id} with SMIRKS: {op_smirks}")
        # print(f"Substrates: {substrates}")
        # print(f"Products: {products}")
        try:
            projected_pdts = project_operator(op_smirks, substrates)
            for proj_smiles in projected_pdts:
                proj_hash = get_molecule_hash(proj_smiles)
                if proj_hash == products_hash:
                    exact_matches.append(op_id)
                    break  # One match is enough
        except Exception as e:
            pass

    return exact_matches

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
        return evodex_data_cache

    # Load the EVODEX data from the JSON file and return it as an object.
    script_dir = os.path.dirname(__file__)
    rel_path = os.path.join('..', 'evodex/data', 'evaluation_operator_data.json')
    json_filepath = os.path.abspath(os.path.join(script_dir, rel_path))
    if os.path.exists(json_filepath):
        with open(json_filepath, 'r') as json_file:
            evodex_data_cache = json.load(json_file)
        return evodex_data_cache
    
    # Index EVODEX data as JSON files, with robust loading for N and C
    e_data = _create_evodex_json('E')
    n_data = {}
    c_data = {}
    try:
        n_data = _create_evodex_json('N')
    except FileNotFoundError as e:
        print(f"Warning: EVODEX-N file not found. N data will be unavailable. {e}")
    try:
        c_data = _create_evodex_json('C')
    except FileNotFoundError as e:
        print(f"Warning: EVODEX-C file not found. C data will be unavailable. {e}")

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
            evodex_dict[source] = {
                "id": evodex_id,
                "smirks": row['smirks']
            }

    with open(json_filepath, 'w') as json_file:
        json.dump(evodex_dict, json_file, indent=4)

    # print(f"EVODEX-{file_suffix} data has been saved to {json_filepath}")
    return evodex_dict

# Example usage:
if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    # methylation_smiles = "CCCO>>CCCOC"
    # is_valid_formula = assign_evodex_F(methylation_smiles)
    # print(f"{methylation_smiles} matches: {is_valid_formula}")

    # matching_operators = match_operators(methylation_smiles, 'E')
    # print(f"Matching operators for {methylation_smiles}: {matching_operators}")

    # exact_matches = find_exact_matching_operators(methylation_smiles, 'E')
    # print(f"Exact matching operators for {methylation_smiles}: {exact_matches}")

    amide_smiles = "CC(=O)NC>>CC(=O)O.NC"
    exact_matches = find_exact_matching_operators(amide_smiles, 'E')
    print(f"Exact matching operators for {amide_smiles}: {exact_matches}")