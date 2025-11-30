#!/usr/bin/env python3

import json
import math
import re
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger
import traceback

RDLogger.DisableLog("rdApp.warning")

from .projectionmap import project_operator_to_mapped_products


def render_reaction_png(reactant, product, out_path, sub_img_size=(300, 300)):
    from rdkit.Chem import Draw

    img = Draw.MolsToGridImage(
        [reactant, product],
        molsPerRow=2,
        subImgSize=sub_img_size,
    )
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


def log_error(
    fh,
    *,
    table_id,
    stage,
    error,
    row_idx=None,
    label=None,
    smiles=None,
    smiles_col=None,
    template=None,
    row_data=None,
    extra=None,
):
    rec = {
        "table_id": table_id,
        "stage": stage,
    }

    # Normalize error information
    if isinstance(error, Exception):
        rec["error_type"] = type(error).__name__
        rec["error_message"] = str(error)
    else:
        rec["error"] = str(error)

    if row_idx is not None:
        rec["row"] = int(row_idx)
    if label is not None:
        rec["label"] = str(label)
    if smiles is not None:
        rec["smiles"] = smiles
        if smiles_col:
            rec[smiles_col] = smiles
    if template is not None:
        rec["template"] = str(template)

    if row_data is not None:
        safe_row = {}
        for k, v in row_data.items():
            safe_row[k] = to_python_value(v)
        rec["row_data"] = safe_row

    if extra:
        rec["extra"] = extra

    if isinstance(error, Exception):
        rec["traceback"] = traceback.format_exc()

    fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    fh.flush()


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
                log_error(
                    err_f,
                    table_id=table_id,
                    stage="input",
                    error=f"missing_csv: {csv_path}",
                )
                continue

            df = pd.read_csv(csv_path)

            smiles_col = pick_smiles_column(df)
            if smiles_col is None:
                log_error(
                    err_f,
                    table_id=table_id,
                    stage="input",
                    error="no_smiles_column_detected",
                )
                continue

            out_jsonl_path = OUT_DIR / f"{table_id}.project.jsonl"
            out_jsonl_path.parent.mkdir(parents=True, exist_ok=True)
            png_table_dir = RENDER_DIR / table_id
            png_table_dir.mkdir(parents=True, exist_ok=True)

            template = meta.get("reaction_template")

            with open(out_jsonl_path, "w", encoding="utf-8") as out_f:
                for idx, row in df.iterrows():
                    smiles_raw = row.get(smiles_col)
                    label = row.get("substrate", f"row_{idx}")
                    row_template = row.get("template") or template

                    if not isinstance(smiles_raw, str) or not smiles_raw.strip():
                        log_error(
                            err_f,
                            table_id=table_id,
                            stage="input",
                            error="missing_or_invalid_substrate_smiles",
                            row_idx=idx,
                            label=label,
                            row_data=row.to_dict(),
                        )
                        continue

                    if not isinstance(row_template, str) or not row_template.strip():
                        log_error(
                            err_f,
                            table_id=table_id,
                            stage="input",
                            error="missing_or_invalid_reaction_template",
                            row_idx=idx,
                            label=label,
                            row_data=row.to_dict(),
                        )
                        continue

                    mol = Chem.MolFromSmiles(smiles_raw)
                    if mol is None:
                        log_error(
                            err_f,
                            table_id=table_id,
                            stage="input",
                            error="invalid_smiles",
                            row_idx=idx,
                            label=label,
                            smiles=smiles_raw,
                            smiles_col=smiles_col,
                            row_data=row.to_dict(),
                        )
                        continue

                    # Projection-based mapping: substrate + operator template
                    try:
                        mapped_substrate_smiles, mapped_product_smiles = project_operator_to_mapped_products(
                            smiles_raw, row_template
                        )
                    except Exception as e:
                        log_error(
                            err_f,
                            table_id=table_id,
                            stage="projection",
                            error=e,
                            row_idx=idx,
                            label=label,
                            smiles=smiles_raw,
                            smiles_col=smiles_col,
                            template=row_template,
                            row_data=row.to_dict(),
                        )
                        continue

                    # For rendering, use the mapped substrate/product so atom maps are visible
                    try:
                        sub_mol = Chem.MolFromSmiles(mapped_substrate_smiles)
                        prod_mol = Chem.MolFromSmiles(mapped_product_smiles)
                    except Exception as e:
                        log_error(
                            err_f,
                            table_id=table_id,
                            stage="render_input",
                            error=f"mol_from_smiles_failed: {e}",
                            row_idx=idx,
                            label=label,
                            smiles=smiles_raw,
                            smiles_col=smiles_col,
                            template=row_template,
                            row_data=row.to_dict(),
                        )
                        continue

                    if sub_mol is None or prod_mol is None:
                        log_error(
                            err_f,
                            table_id=table_id,
                            stage="render_input",
                            error="mol_from_smiles_returned_none",
                            row_idx=idx,
                            label=label,
                            smiles=smiles_raw,
                            smiles_col=smiles_col,
                            template=row_template,
                            row_data=row.to_dict(),
                        )
                        continue

                    slug = slugify(label)
                    png_path = png_table_dir / f"{int(idx):03d}_{slug}.png"
                    try:
                        render_reaction_png(sub_mol, prod_mol, png_path)
                    except Exception as e:
                        log_error(
                            err_f,
                            table_id=table_id,
                            stage="render",
                            error=e,
                            row_idx=idx,
                            label=label,
                            smiles=smiles_raw,
                            smiles_col=smiles_col,
                            template=row_template,
                            row_data=row.to_dict(),
                        )
                        continue

                    # Minimal semantics: raw input substrate + mapped substrate/product + image
                    rec = {
                        "table_id": table_id,
                        "row": int(idx),
                        "label": label,
                        "smiles_input": smiles_raw,
                        "mapped_substrate_smiles": mapped_substrate_smiles,
                        "mapped_product_smiles": mapped_product_smiles,
                        "image": str(png_path.relative_to(BASE_DIR)),
                    }

                    out_f.write(json.dumps(rec, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()