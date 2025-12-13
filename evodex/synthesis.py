import sys
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
import re

# Global cache for EVODEX data (filename -> DataFrame or dict)
evodex_data_cache = {}


def project_reaction_operator(smirks: str, substrate: str):
    """
    Apply a reaction operator (SMIRKS) to a single substrate SMILES.

    Returns a list of unique product SMILES strings.
    """
    rxn = AllChem.ReactionFromSmarts(smirks)
    if not rxn:
        raise ValueError(f"Failed to create reaction from SMIRKS: {smirks}")

    if rxn.GetNumReactantTemplates() == 0 or rxn.GetNumProductTemplates() == 0:
        raise ValueError(f"Reaction has no valid reactants or products: {smirks}")

    substrate_mol = Chem.MolFromSmiles(substrate)
    if not substrate_mol:
        raise ValueError(f"Failed to create molecule from substrate SMILES: {substrate}")

    substrate_mol = Chem.AddHs(substrate_mol)

    products = rxn.RunReactants((substrate_mol,))
    if not products:
        return []

    unique_products = set()
    for product_tuple in products:
        for product in product_tuple:
            if product:
                unique_products.add(Chem.MolToSmiles(product))

    return list(unique_products)


def project_evodex_operator(evodex_id: str, substrate: str):
    """
    Apply an EVODEX operator (by EVODEX id) to a substrate SMILES.
    """
    smirks = _lookup_smirks_by_evodex_id(evodex_id)
    if smirks is None:
        raise ValueError(f"SMIRKS not found for EVODEX ID: {evodex_id}")
    return project_reaction_operator(smirks, substrate)


def project_family_operators(substrate: str, family: str):
    """
    Project all operators for a given EVODEX family on a substrate.

    family examples: "B", "Bm", "C", "Cm", "D", "Dm", "E", "Em"
    """
    ops = _load_evodex_family_operators(family)
    applicable = {}
    for evodex_id, smirks in ops.items():
        try:
            products = project_reaction_operator(smirks, substrate)
            if products:
                applicable[evodex_id] = products
        except Exception:
            pass
    return applicable


def project_synthesis_operators(substrate: str):
    """
    Project synthesis-subset operators on a substrate.

    This keeps the old behavior (subset-only), but uses the new filenames:
    - preferred: EVODEX-D_synthesis_subset.csv
    - fallback:  EVODEX-E_synthesis_subset.csv

    Returns dict: {evodex_id: [product_smiles, ...], ...}
    """
    ops = _load_synthesis_subset_operators()
    applicable = {}
    for evodex_id, smirks in ops.items():
        try:
            products = project_reaction_operator(smirks, substrate)
            if products:
                applicable[evodex_id] = products
        except Exception:
            pass
    return applicable


def _lookup_smirks_by_evodex_id(evodex_id: str):
    """
    Look up SMIRKS by EVODEX ID using new family CSVs:
      evodex/data/EVODEX-<family>.csv
    where <family> is one of: B, C, D, E, Bm, Cm, Dm, Em

    Accepts ids like:
      EVODEX.2-E159
      EVODEX.2-Em159
      EVODEX-E159
      EVODEX-Em159
      EVODEX.1-E159_temp  (suffix ignored; exact id match required)
    """
    family = _infer_family_from_evodex_id(evodex_id)
    if family is None:
        return None

    ops = _load_evodex_family_operators(family)

    if evodex_id in ops:
        return ops[evodex_id]

    # If ids in CSV do not include your suffixes like "_temp", try stripping common suffixes
    base = re.sub(r"(_temp|_draft|_test)$", "", evodex_id)
    return ops.get(base)


def _infer_family_from_evodex_id(evodex_id: str):
    """
    Extract family token after the last "-" in an EVODEX id, allowing an optional "m".

    Examples:
      EVODEX.2-E159      -> "E"
      EVODEX.2-Em159     -> "Em"
      EVODEX.2-D10       -> "D"
      EVODEX.2-Dm10      -> "Dm"
    """
    m = re.search(r"-([BCDE]m?)(\d+)", evodex_id)
    if not m:
        return None
    return m.group(1)


def _load_evodex_family_operators(family: str):
    """
    Load EVODEX-<family>.csv as {id: smirks}.
    """
    filename = f"EVODEX-{family}.csv"
    return _load_evodex_data(filename, key_column="id", value_column="smirks")


def _load_synthesis_subset_operators():
    """
    Load the synthesis subset file if present, preferring D then E.
    """
    for filename in ["EVODEX-D_synthesis_subset.csv", "EVODEX-E_synthesis_subset.csv"]:
        try:
            return _load_evodex_data(filename, key_column="id", value_column="smirks")
        except FileNotFoundError:
            continue
    raise FileNotFoundError(
        "No synthesis subset CSV found. Tried EVODEX-D_synthesis_subset.csv and EVODEX-E_synthesis_subset.csv."
    )


def _load_evodex_data(filename: str, key_column=None, value_column=None):
    """
    Lazy-load EVODEX data from evodex/data/<filename> and cache it.

    If key_column and value_column are provided, returns dict {key: value}.
    Otherwise returns a DataFrame.
    """
    global evodex_data_cache
    if filename not in evodex_data_cache:
        script_dir = os.path.dirname(__file__)
        filepath = os.path.abspath(os.path.join(script_dir, "..", "evodex", "data", filename))
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        df = pd.read_csv(filepath)
        if key_column and value_column:
            evodex_data_cache[filename] = dict(zip(df[key_column], df[value_column]))
        else:
            evodex_data_cache[filename] = df

    return evodex_data_cache[filename]


if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    # Direct projection
    smirks = "[H][C:8]([C:7])([O:9][H])[H:19]>>[C:7][C:8](=[O:9])[H:19]"
    substrate = "CCCO"
    result = project_reaction_operator(smirks, substrate)
    print("Direct projection:", result)

    # Project from EVODEX ID (make sure this id exists in your EVODEX-*.csv)
    evodex_id = "EVODEX.2-E159"
    try:
        result = project_evodex_operator(evodex_id, substrate)
        print("Referenced EVODEX projection:", result)
    except Exception as e:
        print("Referenced EVODEX projection error:", e)

    # Project all operators for a family
    try:
        fam = "E"
        result = project_family_operators(substrate, fam)
        print(f"All EVODEX-{fam} projections (non-empty):", len(result))
    except Exception as e:
        print("Family projection error:", e)

    # Project synthesis subset (if subset CSV exists)
    try:
        result = project_synthesis_operators(substrate)
        print("Synthesis subset projections (non-empty):", len(result))
    except Exception as e:
        print("Synthesis subset projection error:", e)