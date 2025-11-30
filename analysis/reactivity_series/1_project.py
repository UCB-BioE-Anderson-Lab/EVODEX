#!/usr/bin/env python3
"""
Project a set of substrate SMILES through a reaction template and render
mapped substrate and product images.

Inputs
------
- tables_manifest.yaml in the same directory as this script.
  Supported formats:
    1) {"tables": [{"table_id": ..., ...}, ...]}
    2) {"table_id_1": {...}, "table_id_2": {...}, ...}
    3) [{"table_id": ..., ...}, ...]
- CSV files in data/tables/{table_id}.csv
  Each table must contain a substrate SMILES column. The column is detected in
  this order:
    1) "smiles"
    2) "substrate_smiles" or "substrate"
    3) any column name containing "smiles" (case insensitive)
  Each row may optionally contain a "template" column. Otherwise the
  "reaction_template" field from the manifest entry is used.

Outputs
-------
- out/1_project/{table_id}.project.jsonl
    JSON lines file with records:
      {
        "table_id": str,
        "row": int,
        "label": str,
        "smiles_input": str,
        "mapped_substrate_smiles": str,
        "mapped_product_smiles": str,
        "image": str  # path relative to this file's directory
      }

- render/1_project/{table_id}/*.png
    Two panel PNG images showing mapped substrate and product.

- out/_logs/1_project.errors.jsonl
    JSON lines file with error records. Each record contains:
      {
        "table_id": str,
        "stage": str,  # "input", "projection", "render_input", or "render"
        ...
      }
    Depending on the error, additional fields may be present:
      - error_type, error_message or error
      - row (zero based row index)
      - label
      - smiles and original SMILES column name
      - template
      - row_data (sanitized row contents)
      - extra
      - traceback (for exceptions)
"""

import json
import math
import re
import traceback
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger

from .projectionmap import project_operator_to_mapped_products

RDLogger.DisableLog("rdApp.warning")


def render_reaction_png(reactant, product, out_path, sub_img_size=(300, 300)) -> None:
    """
    Render a two panel PNG of reactant and product molecules.

    Parameters
    ----------
    reactant : rdkit.Chem.Mol
        Mapped substrate molecule.
    product : rdkit.Chem.Mol
        Mapped product molecule.
    out_path : pathlib.Path
        Output PNG file path.
    sub_img_size : tuple[int, int], optional
        Per molecule image size in pixels.
    """
    from rdkit.Chem import Draw

    img = Draw.MolsToGridImage(
        [reactant, product],
        molsPerRow=2,
        subImgSize=sub_img_size,
    )
    img.save(str(out_path))


def pick_smiles_column(df: pd.DataFrame):
    """
    Heuristically select the substrate SMILES column from a DataFrame.

    Selection priority
    ------------------
    1) Column named "smiles" (case insensitive)
    2) Column named "substrate_smiles" or "substrate" (case insensitive)
    3) First column whose name contains "smiles" (case insensitive)

    Returns
    -------
    str or None
        The original column name, or None if no candidate is found.
    """
    cols = list(df.columns)

    for col in cols:
        lc = str(col).strip().lower()
        if lc == "smiles":
            return col

    for col in cols:
        lc = str(col).strip().lower()
        if lc in ("substrate_smiles", "substrate"):
            return col

    for col in cols:
        lc = str(col).strip().lower()
        if "smiles" in lc:
            return col

    return None


BASE_DIR = Path(__file__).resolve().parent
TABLES_MANIFEST = BASE_DIR / "tables_manifest.yaml"
TABLES_DIR = BASE_DIR / "data/tables"
OUT_DIR = BASE_DIR / "out/1_project"
RENDER_DIR = BASE_DIR / "render/1_project"
LOG_DIR = BASE_DIR / "out/_logs"

# Optional manual filter for debugging specific tables.
ONLY_TABLES = []  # for example: ["trifluoroacetanilide_methanolysis"]


def load_manifest(path: Path) -> dict:
    """
    Load the tables manifest and return a mapping from table_id to metadata.

    Supported formats
    -----------------
    1) {"tables": [{"table_id": ..., ...}, ...]}
    2) {"table_id_1": {...}, "table_id_2": {...}, ...}
       (table_id is taken from the key if missing inside the record)
    3) [{"table_id": ..., ...}, ...]

    Parameters
    ----------
    path : pathlib.Path
        Path to the YAML manifest.

    Returns
    -------
    dict[str, dict]
        Mapping from table_id to manifest record.
    """
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
    """
    Convert an arbitrary label to a filesystem friendly slug.

    Rules
    -----
    - Strip leading and trailing whitespace
    - Lowercase
    - Replace internal whitespace runs with a single underscore
    - Remove all characters except a-z, 0-9 and underscore
    - Use "entry" if the result is empty

    Parameters
    ----------
    label : str
        Input label.

    Returns
    -------
    str
        Slugified label.
    """
    s = str(label).strip().lower()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^a-z0-9_]+", "", s)
    if not s:
        s = "entry"
    return s


def to_python_value(v):
    """
    Convert a value to a JSON safe plain Python type.

    Behavior
    --------
    - If v has an .item attribute (for example a NumPy scalar), use v.item().
    - On failure to access .item, fall back to str(v).
    - Replace NaN or infinite floats with None.

    Parameters
    ----------
    v : any
        Input value.

    Returns
    -------
    any
        Converted value suitable for JSON encoding.
    """
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
) -> None:
    """
    Write a single error record as JSON to an open file handle.

    Parameters
    ----------
    fh : io.TextIOBase
        Open file handle for the error log.
    table_id : str
        Identifier of the table being processed.
    stage : str
        Pipeline stage where the error occurred.
        For example: "input", "projection", "render_input", "render".
    error : Exception or str
        Exception instance or error description.
    row_idx : int, optional
        Zero based index of the failing row in the source table.
    label : str, optional
        Human readable label for the entry.
    smiles : str, optional
        Substrate SMILES string that led to the error.
    smiles_col : str, optional
        Name of the SMILES column.
    template : str, optional
        Reaction template used for projection.
    row_data : dict, optional
        Raw row contents from the source CSV.
    extra : dict, optional
        Any additional diagnostic information.
    """
    rec = {
        "table_id": table_id,
        "stage": stage,
    }

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


def main() -> None:
    """
    Run the projection and rendering pipeline for all tables in the manifest.

    For each table:
      - Load the corresponding CSV file.
      - Detect the substrate SMILES column.
      - Resolve the reaction template per row or from the manifest.
      - Project mapped substrate and product SMILES.
      - Render the mapped molecules as a PNG image.
      - Write JSON records describing the successful projections.

    All failures are recorded in out/_logs/1_project.errors.jsonl.
    """
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

                    try:
                        mapped_substrate_smiles, mapped_product_smiles = (
                            project_operator_to_mapped_products(
                                smiles_raw,
                                row_template,
                            )
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