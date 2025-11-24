#!/usr/bin/env python3

import json
import math
import re
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem


from rdkit.Chem import rdChemReactions

def clear_isotopes(mol):
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return mol

def sanitize_lenient(mol):
    try:
        Chem.SanitizeMol(mol)
        return mol, None
    except Exception as e:
        return None, str(e)

def fix_overvalent_nitrogens(mol):
    try:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1 and atom.GetTotalValence() > 4:
                atom.SetFormalCharge(0)
        Chem.SanitizeMol(mol)
        return mol, None
    except Exception as e:
        return None, str(e)

def to_explicit_h_smiles(mol):
    mol2 = Chem.AddHs(mol)
    return Chem.MolToSmiles(mol2, kekuleSmiles=False)

def apply_template_first(template, reactant_mol):
    try:
        rxn = rdChemReactions.ReactionFromSmarts(template)
        if rxn is None:
            return None, "ReactionFromSmarts returned None"
        outcomes = rxn.RunReactants((reactant_mol,))
        if not outcomes:
            return None, "no_products"
        prod = Chem.Mol(outcomes[0][0])
        prod, err = sanitize_lenient(prod)
        if prod is None:
            repaired, rerr = fix_overvalent_nitrogens(Chem.Mol(outcomes[0][0]))
            if repaired is None:
                return None, f"product_{err}; repair={rerr}"
            prod = repaired
        prod = clear_isotopes(prod)
        return prod, None
    except Exception as e:
        return None, f"reaction_failed: {e}"

def render_reaction_png(reactant, product, out_path, sub_img_size=(300,300)):
    from rdkit.Chem import Draw
    img = Draw.MolsToGridImage([reactant, product], molsPerRow=2, subImgSize=sub_img_size)
    img.save(str(out_path))

def pick_smiles_column(df):
    cols = list(df.columns)
    for c in cols:
        lc = str(c).strip().lower()
        if lc == "smiles":
            return c
    for c in cols:
        lc = str(c).strip().lower()
        if lc in ("substrate_smiles", "substrate"):
            return c
    for c in cols:
        lc = str(c).strip().lower()
        if "smiles" in lc:
            return c
    return None


BASE_DIR = Path(__file__).resolve().parent
TABLES_MANIFEST = BASE_DIR / "tables_manifest.yaml"
TABLES_DIR = BASE_DIR / "data/tables"
OUT_DIR = BASE_DIR / "out/1_project"
RENDER_DIR = BASE_DIR / "render/1_project"
LOG_DIR = BASE_DIR / "out/_logs"

ONLY_TABLES = []  # e.g. ["trifluoroacetanilide_methanolysis"]


def load_manifest(path: str):
    with open(path, "r", encoding="utf-8") as f:
        obj = yaml.safe_load(f)

    entries = []

    if isinstance(obj, dict) and "tables" in obj:
        for rec in obj["tables"]:
            entries.append(rec)
    elif isinstance(obj, dict):
        for table_id, rec in obj.items():
            if isinstance(rec, dict):
                rec = dict(rec)
                rec.setdefault("table_id", table_id)
                entries.append(rec)
    elif isinstance(obj, list):
        entries = obj
    else:
        raise ValueError("Unsupported manifest format")

    by_id = {}
    for rec in entries:
        tid = rec.get("table_id")
        if not tid:
            continue
        by_id[tid] = rec
    return by_id


def slugify(label: str) -> str:
    s = str(label).strip().lower()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^a-z0-9_]+", "", s)
    if not s:
        s = "entry"
    return s


def to_python_value(v):
    try:
        if hasattr(v, "item"):
            v = v.item()
    except Exception:
        v = str(v)
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return None
    return v


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RENDER_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(TABLES_MANIFEST)
    errors_path = LOG_DIR / "1_project.errors.jsonl"

    with open(errors_path, "w", encoding="utf-8") as err_f:
        for table_id, meta in sorted(manifest.items()):
            if ONLY_TABLES and table_id not in ONLY_TABLES:
                continue

            csv_path = TABLES_DIR / f"{table_id}.csv"
            if not csv_path.exists():
                err_rec = {
                    "table_id": table_id,
                    "stage": "input",
                    "error": f"missing_csv: {csv_path}",
                }
                err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                continue

            df = pd.read_csv(csv_path)

            smiles_col = pick_smiles_column(df)

            value_cols = [
                c for c in df.columns
                if c not in (smiles_col, "template", "label", "smiles")
            ]

            out_jsonl_path = OUT_DIR / f"{table_id}.project.jsonl"
            out_jsonl_path.parent.mkdir(parents=True, exist_ok=True)
            png_table_dir = RENDER_DIR / table_id
            png_table_dir.mkdir(parents=True, exist_ok=True)

            template = meta.get("reaction_template")
            with open(out_jsonl_path, "w", encoding="utf-8") as out_f:
                for idx, row in df.iterrows():
                    smiles_raw = row.get(smiles_col)
                    label = row.get("label", f"row_{idx}")
                    row_template = row.get("template", template)

                    if not isinstance(smiles_raw, str) or not smiles_raw.strip():
                        err_rec = {
                            "table_id": table_id,
                            "row": int(idx),
                            "label": label,
                            "stage": "input",
                            "error": "missing_or_invalid_substrate_smiles",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    mol = Chem.MolFromSmiles(smiles_raw)
                    if mol is None:
                        err_rec = {
                            "table_id": table_id,
                            "row": int(idx),
                            "label": label,
                            smiles_col: smiles_raw,
                            "stage": "input",
                            "error": "invalid_smiles",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    prod, perr = apply_template_first(row_template, mol)
                    if prod is None:
                        err_rec = {
                            "table_id": table_id,
                            "row": int(idx),
                            "label": label,
                            smiles_col: smiles_raw,
                            "stage": "reaction",
                            "error": perr or "no_products",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    substrate_canon = Chem.MolToSmiles(mol, kekuleSmiles=False)
                    prod_smiles = Chem.MolToSmiles(prod, kekuleSmiles=False)
                    prod_smiles_h = to_explicit_h_smiles(prod)

                    slug = slugify(label)
                    png_path = png_table_dir / f"{int(idx):03d}_{slug}.png"
                    try:
                        render_reaction_png(mol, prod, png_path)
                    except Exception as e:
                        err_rec = {
                            "table_id": table_id,
                            "row": int(idx),
                            "label": label,
                            smiles_col: smiles_raw,
                            "stage": "render",
                            "error": f"{e}",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    row_data = {}
                    for col in value_cols:
                        row_data[col] = to_python_value(row.get(col))

                    rec = {
                        "table_id": table_id,
                        "row": int(idx),
                        "label": label,
                        f"{smiles_col}_input": smiles_raw,
                        smiles_col: substrate_canon,
                        "product_smiles": prod_smiles,
                        "product_smiles_explicit_h": prod_smiles_h,
                        "image": str(png_path.relative_to(BASE_DIR)),
                        "reaction_template": row_template,
                        "reactivity_type": meta.get("reactivity_type"),
                        "units": meta.get("units"),
                        "context": meta.get("context"),
                        "source_ref": meta.get("source_ref"),
                        "data": row_data,
                    }

                    out_f.write(json.dumps(rec, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()