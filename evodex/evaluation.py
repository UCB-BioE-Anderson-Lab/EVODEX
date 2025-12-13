import os
import sys
import json
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import AllChem
from itertools import combinations, permutations

from evodex.formula import calculate_formula_diff
from evodex.utils import get_molecule_hash_from_mol
from evodex.projection import project_operator


# Initialize caches
evodex_f_cache = None
evodex_data_cache = None


# Back-compat helper used throughout this file
def get_molecule_hash(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for hashing: {smiles}")
    return get_molecule_hash_from_mol(mol)


def _add_hydrogens(rxn_smiles: str) -> str:
    """
    Add explicit Hs to both sides of a reaction SMILES (substrates>>products).
    """
    substrate, product = rxn_smiles.split(">>")
    substrate_mol = Chem.MolFromSmiles(substrate)
    product_mol = Chem.MolFromSmiles(product)
    if substrate_mol is None or product_mol is None:
        raise ValueError(f"Invalid reaction SMILES: {rxn_smiles}")

    substrate_mol = Chem.AddHs(substrate_mol)
    product_mol = Chem.AddHs(product_mol)

    substrate_smiles = Chem.MolToSmiles(substrate_mol)
    product_smiles = Chem.MolToSmiles(product_mol)
    return f"{substrate_smiles}>>{product_smiles}"


def operator_matches_reaction(operator_smirks: str, reaction_smiles: str) -> bool:
    """
    True if operator_smirks can produce reaction_smiles products from its substrates.
    reaction_smiles must be substrates>>products in SMILES.
    """
    if not isinstance(operator_smirks, str) or not operator_smirks.strip():
        raise ValueError(f"[operator_matches_reaction] Invalid operator_smirks: {operator_smirks}")
    if not isinstance(reaction_smiles, str) or not reaction_smiles.strip():
        raise ValueError(f"[operator_matches_reaction] Invalid reaction_smiles: {reaction_smiles}")

    op_sub_parts = operator_smirks.split(">>")[0].split(".")
    rxn_sub_parts = reaction_smiles.split(">>")[0].split(".")
    if len(rxn_sub_parts) != len(op_sub_parts):
        return False

    op_prod_parts = operator_smirks.split(">>")[1].split(".")
    rxn_prod_parts = reaction_smiles.split(">>")[1].split(".")
    if len(rxn_prod_parts) != len(op_prod_parts):
        return False

    op_sub_smarts = [Chem.MolFromSmarts(m) for m in op_sub_parts]
    if any(m is None for m in op_sub_smarts):
        return False

    rxn_sub_mols = []
    for m in rxn_sub_parts:
        mol = Chem.MolFromSmiles(m)
        if mol is None:
            return False
        rxn_sub_mols.append(Chem.AddHs(mol))

    found_match = False
    for op_perm in permutations(op_sub_smarts):
        used = set()
        ok = True
        for op_mol in op_perm:
            matched_one = False
            for i, rxn_mol in enumerate(rxn_sub_mols):
                if i in used:
                    continue
                if rxn_mol.HasSubstructMatch(op_mol):
                    used.add(i)
                    matched_one = True
                    break
            if not matched_one:
                ok = False
                break
        if ok:
            found_match = True
            break

    if not found_match:
        return False

    try:
        reaction_smiles_with_h = _add_hydrogens(reaction_smiles)

        rxn = AllChem.ReactionFromSmarts(reaction_smiles_with_h, useSmiles=True)
        rxn.Initialize()

        # Strip atom maps from templates
        for i in range(rxn.GetNumReactantTemplates()):
            for atom in rxn.GetReactantTemplate(i).GetAtoms():
                atom.SetAtomMapNum(0)
        for i in range(rxn.GetNumProductTemplates()):
            for atom in rxn.GetProductTemplate(i).GetAtoms():
                atom.SetAtomMapNum(0)

        substrates = ".".join(
            Chem.MolToSmiles(rxn.GetReactantTemplate(i)) for i in range(rxn.GetNumReactantTemplates())
        )
        expected_products = ".".join(
            Chem.MolToSmiles(rxn.GetProductTemplate(i)) for i in range(rxn.GetNumProductTemplates())
        )

        expected_hashes = {get_molecule_hash(p) for p in expected_products.split(".")}

        projected_product_sets = project_operator(operator_smirks, substrates)
        for product_string in projected_product_sets:
            projected_hashes = {get_molecule_hash(p) for p in product_string.split(".")}
            if projected_hashes == expected_hashes:
                return True

        return False
    except Exception:
        return False


def assign_evodex_F(rxn_smiles: str):
    """
    Assign EVODEX-F id(s) to a reaction SMILES using formula-diff matching.

    Returns:
      - list[str] of matching F ids (may be empty)
    """
    smirks_with_h = _add_hydrogens(rxn_smiles)
    formula_diff = calculate_formula_diff(smirks_with_h)
    evodex_f = _load_evodex_f()
    return evodex_f.get(frozenset(formula_diff.items()), [])


def _load_evodex_f():
    """
    Load EVODEX-F cache from EVODEX-F.csv.

    Cache format:
      frozenset(formula_diff.items()) -> [F_id, F_id, ...]
    """
    global evodex_f_cache
    if evodex_f_cache is not None:
        return evodex_f_cache

    evodex_f_cache = {}
    script_dir = os.path.dirname(__file__)
    filepath = os.path.abspath(os.path.join(script_dir, "..", "evodex", "data", "EVODEX-F.csv"))
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    evodex_f_df = pd.read_csv(filepath)
    for _, row in evodex_f_df.iterrows():
        formula_diff = eval(row["formula_diff"])
        f_id = row["id"]
        key = frozenset(formula_diff.items())
        evodex_f_cache.setdefault(key, []).append(f_id)

    return evodex_f_cache


def _parse_sources(sources: str):
    sources = sources.replace('"', "")
    return sources.split(",")


def match_operators(rxn_smiles: str, evodex_type: str = "E"):
    """
    Given a reaction SMILES (possibly multi-reactant/product), try all subset pairings
    and return all operator IDs that exactly produce the paired product(s).

    evodex_type is now one of:
      B, Bm, C, Cm, D, Dm, E, Em
    """
    valid_operators = []

    if ">>" not in rxn_smiles:
        return valid_operators

    try:
        substrates, products = rxn_smiles.split(">>")
        substrate_list = substrates.split(".")
        product_list = products.split(".")

        substrate_indices = list(range(len(substrate_list)))
        product_indices = list(range(len(product_list)))

        all_pairings = set()
        for i in range(1, len(substrate_indices) + 1):
            for j in range(1, len(product_indices) + 1):
                for reactant_combo in combinations(substrate_indices, i):
                    for product_combo in combinations(product_indices, j):
                        all_pairings.add((frozenset(reactant_combo), frozenset(product_combo)))

        for reactant_idx, product_idx in all_pairings:
            reactant_smiles = ".".join(substrate_list[i] for i in sorted(reactant_idx))
            product_smiles = ".".join(product_list[i] for i in sorted(product_idx))
            pairing_smiles = f"{reactant_smiles}>>{product_smiles}"
            valid_operators.extend(_match_operator(pairing_smiles, evodex_type))

    except Exception as e:
        print(f"Error processing SMILES {rxn_smiles}: {e}")

    return valid_operators


def _match_operator(rxn_smiles: str, evodex_type: str = "E"):
    """
    For a single paired reaction SMILES (reactants>>products), do:
      formula-diff -> candidate F ids -> candidate operators in that family -> exact projection match

    Returns: list[str] operator IDs
    """
    smiles_with_h = _add_hydrogens(rxn_smiles)
    formula_diff = calculate_formula_diff(smiles_with_h)

    evodex_f = _load_evodex_f()
    f_id_list = evodex_f.get(frozenset(formula_diff.items()), [])
    if not f_id_list:
        return []

    f_id = f_id_list[0]  # keep old behavior: pick first

    evodex_data = _load_evodex_data()
    if f_id not in evodex_data:
        return []

    potential_operators = evodex_data[f_id].get(evodex_type, [])
    if not potential_operators:
        return []

    sub_smiles, pdt_smiles = rxn_smiles.split(">>")
    pdt_hash = get_molecule_hash(pdt_smiles)

    valid = []
    for operator in potential_operators:
        try:
            projected_pdts = project_operator(operator["smirks"], sub_smiles)
            for proj_smiles in projected_pdts:
                if get_molecule_hash(proj_smiles) == pdt_hash:
                    valid.append(operator["id"])
                    break
        except Exception:
            pass

    return valid


def find_exact_matching_operators(p_smiles: str, evodex_type: str = "E"):
    """
    Find operators that exactly match the given P reaction SMILES (substrates>>products).

    evodex_type is now one of:
      B, Bm, C, Cm, D, Dm, E, Em
    """
    f_id_list = assign_evodex_F(p_smiles)
    if not f_id_list:
        return []

    f_id = f_id_list[0]

    evodex_data = _load_evodex_data()
    if f_id not in evodex_data:
        return []

    candidate_operators = evodex_data[f_id].get(evodex_type, [])
    if not candidate_operators:
        return []

    substrates, products = p_smiles.split(">>")
    products_hash = get_molecule_hash(products)

    exact = []
    for operator in candidate_operators:
        op_id = operator["id"]
        op_smirks = operator.get("smirks")
        if not op_smirks:
            continue
        try:
            projected_pdts = project_operator(op_smirks, substrates)
            for proj_smiles in projected_pdts:
                if get_molecule_hash(proj_smiles) == products_hash:
                    exact.append(op_id)
                    break
        except Exception:
            pass

    return exact


def _load_evodex_data():
    """
    Build or load a JSON index:
      F_id -> { family_token -> [ {id, smirks}, ... ] }

    Families supported:
      B, Bm, C, Cm, D, Dm, E, Em
    """
    global evodex_data_cache
    if evodex_data_cache is not None:
        return evodex_data_cache

    script_dir = os.path.dirname(__file__)
    json_filepath = os.path.abspath(os.path.join(script_dir, "..", "evodex", "data", "evaluation_operator_data.json"))

    if os.path.exists(json_filepath):
        with open(json_filepath, "r") as f:
            evodex_data_cache = json.load(f)
        return evodex_data_cache

    # Build per-family mapping: P_id -> {id, smirks}
    family_tokens = ["B", "Bm", "C", "Cm", "D", "Dm", "E", "Em"]
    family_data = {}
    for fam in family_tokens:
        try:
            family_data[fam] = _create_evodex_json(fam)
        except FileNotFoundError as e:
            print(f"Warning: EVODEX-{fam}.csv not found; family unavailable. {e}")
            family_data[fam] = {}

    evodex_data_cache = {}

    # Load EVODEX-F to map formula diff to P sources
    f_csv = os.path.abspath(os.path.join(script_dir, "..", "evodex", "data", "EVODEX-F.csv"))
    evodex_f_df = pd.read_csv(f_csv)

    for _, row in evodex_f_df.iterrows():
        f_id = row["id"]
        p_ids = _parse_sources(row["sources"])

        per_f = {fam: [] for fam in family_tokens}

        for p_id in p_ids:
            for fam in family_tokens:
                if p_id in family_data[fam]:
                    op = family_data[fam][p_id]
                    if op not in per_f[fam]:
                        per_f[fam].append(op)

        evodex_data_cache[f_id] = per_f

    with open(json_filepath, "w") as f:
        json.dump(evodex_data_cache, f, indent=4)

    return evodex_data_cache


def _create_evodex_json(family_token: str):
    """
    Create a dict mapping each P source id -> {id, smirks} for a given operator family.

    Reads:
      evodex/data/EVODEX-<family_token>.csv

    Writes (optional cache):
      evodex/data/evodex_<family_token.lower()>_data.json
    """
    script_dir = os.path.dirname(__file__)
    csv_filepath = os.path.abspath(os.path.join(script_dir, "..", "evodex", "data", f"EVODEX-{family_token}.csv"))
    json_filepath = os.path.abspath(os.path.join(script_dir, "..", "evodex", "data", f"evodex_{family_token.lower()}_data.json"))

    if not os.path.exists(csv_filepath):
        raise FileNotFoundError(f"File not found: {csv_filepath}")

    df = pd.read_csv(csv_filepath)

    evodex_dict = {}
    for _, row in df.iterrows():
        evodex_id = row["id"]
        sources = _parse_sources(row["sources"])
        smirks = row["smirks"]
        for source in sources:
            evodex_dict[source] = {"id": evodex_id, "smirks": smirks}

    with open(json_filepath, "w") as f:
        json.dump(evodex_dict, f, indent=4)

    return evodex_dict


if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    # smiles = "CC(=O)NC>>CC(=O)O.NC" # amide hydrolysis
    # smiles = "CCCO>>CCCOC"  # methylation
    smiles = "CCCO>>CCC=O" # oxidation

    exact_matches = find_exact_matching_operators(smiles, "E")
    print(f"Exact matching operators for {smiles}: {exact_matches}")