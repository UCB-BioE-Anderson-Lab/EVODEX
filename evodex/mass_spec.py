import os
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

import pandas as pd

from evodex.synthesis import project_evodex_operator
from evodex.evaluation import _load_evodex_data, _parse_sources
from evodex.utils import get_molecule_hash_from_mol


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
    raise ValueError(f"Invalid SMILES string: {smiles}")


def _evodex_data_dir():
    """
    Return absolute path to evodex/data directory.
    mass_spec.py lives in EVODEX/evodex/, so data is EVODEX/evodex/data/.
    """
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))


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

        filepath = os.path.join(_evodex_data_dir(), "EVODEX-M_mass_spec_subset.csv")
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        df = pd.read_csv(filepath)

        for _, row in df.iterrows():
            evodex_m_cache.append(
                {
                    "id": row["id"],
                    "mass": float(row["mass"]),
                    "sources": _parse_sources(row["sources"]),
                }
            )

    return evodex_m_cache


def _load_evodex_m_to_f():
    """
    Load or create EVODEX-M -> EVODEX-F mapping.

    EVODEX-M has sources that are EVODEX-P IDs.
    EVODEX-F.csv has sources that are EVODEX-P IDs.
    So: M -> (P sources) -> (F IDs).

    Returns:
    dict: { evodex_m_id: [evodex_f_id, ...] }
    """
    global evodex_m_to_f_cache
    if evodex_m_to_f_cache is not None:
        return evodex_m_to_f_cache

    mapping_path = os.path.join(_evodex_data_dir(), "evodex_m_to_F_mapping.csv")

    if os.path.exists(mapping_path):
        evodex_m_to_f_cache = {}
        df = pd.read_csv(mapping_path)
        for _, row in df.iterrows():
            evodex_m_to_f_cache.setdefault(row["evodex_m"], []).append(row["evodex_f"])
        return evodex_m_to_f_cache

    # Load EVODEX-M
    evodex_m_cache_local = _load_evodex_m()

    # Build P -> F from EVODEX-F.csv
    f_path = os.path.join(_evodex_data_dir(), "EVODEX-F.csv")
    if not os.path.exists(f_path):
        raise FileNotFoundError(f"File not found: {f_path}")

    evodex_f_df = pd.read_csv(f_path)
    p_to_f_map = {}
    for _, row in evodex_f_df.iterrows():
        f_id = row["id"]
        p_ids = _parse_sources(row["sources"])
        for p_id in p_ids:
            p_to_f_map.setdefault(p_id, []).append(f_id)

    # M -> F via shared P sources
    evodex_m_to_f_cache = {}
    for entry in evodex_m_cache_local:
        evodex_m_id = entry["id"]
        for p_id in entry["sources"]:
            if p_id in p_to_f_map:
                evodex_m_to_f_cache.setdefault(evodex_m_id, []).extend(p_to_f_map[p_id])

    # Save mapping for reuse
    with open(mapping_path, "w") as f:
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
    return [
        {"id": entry["id"], "mass": entry["mass"]}
        for entry in evodex_m
        if abs(entry["mass"] - mass_diff) <= precision
    ]


def get_reaction_operators(mass_diff, precision=0.01):
    """
    Retrieve operators (by family) that could explain the mass difference.

    Returns:
    tuple:
        matching_operators: dict {op_type: [operator_dict, ...]} for all families present
        matching_evodex_m: list of matching EVODEX-M rows
        f_ids: sorted list of EVODEX-F IDs implicated
    """
    matching_evodex_m = find_evodex_m(mass_diff, precision)
    if not matching_evodex_m:
        return {}, matching_evodex_m, []

    evodex_m_to_f = _load_evodex_m_to_f()
    evodex_data = _load_evodex_data()

    matching_operators = {}  # dynamic: B/C/D/E/etc
    f_ids = set()

    # Debug: track which families appear
    seen_families = set()

    for entry in matching_evodex_m:
        evodex_m_id = entry["id"]
        for f_id in evodex_m_to_f.get(evodex_m_id, []):
            f_ids.add(f_id)

            if f_id not in evodex_data:
                continue

            blocks = evodex_data[f_id]
            if not isinstance(blocks, dict):
                continue

            for op_type, ops in blocks.items():
                # ignore metadata keys or unexpected payloads
                if not isinstance(ops, list):
                    continue

                seen_families.add(op_type)
                matching_operators.setdefault(op_type, []).extend(ops)

    if seen_families:
        print(f"[mass_spec] operator families seen: {sorted(seen_families)}")

    return matching_operators, matching_evodex_m, sorted(f_ids)


def predict_products(smiles, mass_diff, precision=0.01, evodex_type="E"):
    """
    Project EVODEX operators of a selected family consistent with the EVODEX-M onto a given substrate.

    Parameters:
    smiles (str): substrate SMILES.
    mass_diff (float): target delta mass.
    precision (float): matching tolerance.
    evodex_type (str): operator family to project (e.g. 'B','C','D','E').

    Returns:
    dict: { product_hash: {smiles: str, projections: {(f_ids_tuple, m_id): [op_id,...]}}}
    """
    matching_operators, matching_evodex_m, f_ids = get_reaction_operators(mass_diff, precision)
    if not matching_evodex_m:
        return {}

    if evodex_type not in matching_operators:
        available = sorted(matching_operators.keys())
        print(f"[PRED] requested evodex_type={evodex_type} not found. Available: {available}")
        return {}

    evodex_ops = matching_operators[evodex_type]
    print(f"[PRED] substrate={smiles}")
    print(f"[PRED] evodex_type={evodex_type}, ops={len(evodex_ops)}")

    results = {}
    for operator in evodex_ops:
        try:
            m_id = matching_evodex_m[0]["id"]
            operator_id = operator["id"]

            projected_pdts = project_evodex_operator(operator_id, smiles)

            for proj_smiles in projected_pdts:
                mol = Chem.MolFromSmiles(proj_smiles)
                if mol is None:
                    continue

                proj_hash = get_molecule_hash_from_mol(mol)

                if proj_hash not in results:
                    results[proj_hash] = {"smiles": proj_smiles, "projections": {}}

                formula_mass_key = (tuple(f_ids), m_id)
                results[proj_hash]["projections"].setdefault(formula_mass_key, [])

                if operator_id not in results[proj_hash]["projections"][formula_mass_key]:
                    results[proj_hash]["projections"][formula_mass_key].append(operator_id)

        except Exception as e:
            print(f"{operator.get('id', '<noid>')} errored: {str(e)}")

    return results


if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    substrate = "CCCO"
    mass_diff = 14.016
    precision = 0.01

    evodex_m = find_evodex_m(mass_diff, precision)
    print(f"Found matching {mass_diff}: {evodex_m}")

    matching_operators, _, f_ids = get_reaction_operators(mass_diff, precision)
    print(f"F IDs: {f_ids}")

    # Print operator counts by family (B/C/D/E/etc)
    fam_counts = {k: len(v) for k, v in matching_operators.items()}
    print(f"Operator counts by family: {fam_counts}")

    # Example projection: change this to 'B', 'C', 'D', or 'E'
    evodex_type = "E"

    results = predict_products(substrate, mass_diff, precision, evodex_type)
    for _, details in results.items():
        print(f"Product: {details['smiles']}")
        for (f_ids_key, m_id), operators in details["projections"].items():
            print(f"  EVODEX-F IDs: {f_ids_key}, EVODEX-M ID: {m_id}, Operators: {operators}")