from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from evodex.utils import get_molecule_hash
from evodex.formula import calculate_formula_diff
from typing import List, Dict
import itertools

def project_operator(operator_smirks, substrates):
    """
    Apply a reaction operator (SMIRKS) to an n-substrate set and return a list of product SMILES sets.
    """

    # Check type and format of operator_smirks
    if not isinstance(operator_smirks, str) or '>>' not in operator_smirks:
        print(f"[project_operator] Invalid operator_smirks: {operator_smirks}")
        return []
    # Build reaction from SMIRKS
    try:
        substrate_smirks, product_smirks = operator_smirks.split('>>')
        reactant_templates = substrate_smirks.split('.')
        num_operator_substrates = len(reactant_templates)
        rxn = AllChem.ChemicalReaction()
        for r in reactant_templates:
            mol = Chem.MolFromSmarts(r)
            if mol is None:
                return []
            rxn.AddReactantTemplate(mol)

        product_template = Chem.MolFromSmarts(product_smirks)
        if product_template is None:
            return []
        rxn.AddProductTemplate(product_template)
    except Exception as e:
        print(f"Error parsing operator SMIRKS: {e}")
        return []

    # Convert input substrates into molecules
    smiles_list = substrates.split('.')
    substrate_mols = []
    for smi in smiles_list:
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        if mol is None:
            return []
        substrate_mols.append(mol)

    num_input_substrates = len(substrate_mols)

    # Align input substrates with operator requirements
    if num_operator_substrates == 1 and num_input_substrates == 1:
        permutations_to_try = [tuple(substrate_mols)]
    elif num_operator_substrates == num_input_substrates:
        permutations_to_try = list(itertools.permutations(substrate_mols, num_operator_substrates))
    elif num_operator_substrates > 1 and num_input_substrates == 1:
        repeated = [substrate_mols[0]] * num_operator_substrates
        permutations_to_try = [tuple(repeated)]
    else:
        return []

    all_products = []
    for perm_group in permutations_to_try:
        try:
            result = rxn.RunReactants(tuple(perm_group))
            if result:
                all_products.extend(result)
        except Exception as e:
            pass

    unique_products = set()
    for product_tuple in all_products:
        product_smiles = [Chem.MolToSmiles(p) for p in product_tuple if p]
        joined = '.'.join(sorted(product_smiles))
        unique_products.add(joined)

    return list(unique_products)

def match_projection(ero_smirks, substrate, expected_product):
    """
    Project the operator and compare against expected product. Return True if matched.
    """
    try:
        projected = project_operator(ero_smirks, substrate)
        target_hash = get_molecule_hash(expected_product)
        return any(get_molecule_hash(p) == target_hash for p in projected)
    except Exception:
        return False

def add_explicit_hydrogens(smirks):
    """
    Add hydrogens to both sides of a SMIRKS string.
    """
    try:
        substrate, product = smirks.split('>>')
        sub_mol = Chem.AddHs(Chem.MolFromSmiles(substrate))
        prod_mol = Chem.AddHs(Chem.MolFromSmiles(product))
        return f"{Chem.MolToSmiles(sub_mol)}>>{Chem.MolToSmiles(prod_mol)}"
    except Exception as e:
        raise ValueError(f"Could not add hydrogens to SMIRKS: {smirks}") from e

def compute_formula_difference(smirks):
    """
    Compute formula difference after hydrogen normalization.
    """
    smirks_h = add_explicit_hydrogens(smirks)
    return calculate_formula_diff(smirks_h)

def find_matching_eros(evop_smirks: str, candidate_eros: List[Dict]) -> List[Dict]:
    """
    Given a reaction SMIRKS (typically from an EVODEX-P reaction) and a list of candidate
    reaction operators (EROs), identify which operators successfully project the correct
    product from the given substrate.

    Parameters:
    evop_smirks (str): The SMIRKS string representing the input reaction to explain.
    candidate_eros (List[Dict]): A list of EROs where each ERO is a dictionary that includes
        at least 'id' and 'smirks'. The EROs are assumed to have been filtered by formula
        difference (EVODEX-F), though this is not strictly required.

    Returns:
    List[Dict]: A list of matching EROs. Each entry is a dictionary containing:
        - 'id': The ID of the matching ERO.
        - 'smirks': The SMIRKS string of the matching ERO.
        - 'matched_smiles': The projected product that matches the EVOP.
        - 'label_map': A dictionary mapping atom indices from the EVOP to the matching ERO.
                       (Currently placeholder; should be implemented when mapping is available.)
        - 'ero_hash': A unique determinate hash for the reaction operator
    """
    try:
        substrate, product = evop_smirks.split('>>')
        expected_hash = get_molecule_hash(product)
    except Exception as e:
        raise ValueError(f"Invalid EVOP SMIRKS: {evop_smirks}") from e

    matched = []
    for ero in candidate_eros:
        try:
            projected = project_operator(ero['smirks'], substrate)
            for p in projected:
                if get_molecule_hash(p) == expected_hash:
                    matched.append({
                        'id': ero['id'],
                        'smirks': ero['smirks'],
                        'matched_smiles': p,
                        'label_map': {},  # TODO: populate with actual atom mapping
                        'ero_hash': ero['ero_hash']
                    })
                    break
        except Exception:
            continue
    return matched