

#!/usr/bin/env python3
"""
Iterate over reactivity tables listed in a YAML manifest, apply an RDKit
reaction SMARTS template to each substrate to generate products, and render an
image of the resulting reaction instance (reactant -> product).

Adds robust error tracking so you can see exactly which CSV and row failed
(e.g., explicit valence or kekulization issues). Per-table `errors.csv` files are written.

Requires: rdkit, pandas, pyyaml
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import yaml

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem import rdDepictor

RDLogger.DisableLog("rdApp.error")  # silence C++ noise; we log context ourselves

# -----------------------------
# Helpers
# -----------------------------

def slugify(text: str, max_len: int = 40) -> str:
    text = re.sub(r"[^A-Za-z0-9\-_. ]+", "", str(text)).strip().lower()
    text = re.sub(r"[\s/]+", "_", text)
    return (text[:max_len] or "row").strip("._")


def load_manifest(path: Path) -> List[Dict]:
    with open(path, "r") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, list):
        raise ValueError("Manifest YAML must be a list of table entries.")
    required = {"table_id", "reaction_template"}
    for entry in data:
        missing = required - set(entry)
        if missing:
            raise ValueError(f"Manifest entry missing required keys: {missing}")
    return data


_PREFERRED_SMILES_COLS = [
    "smiles",
    "substrate_smiles",
    "reactant_smiles",
    "substrate",
    "reactant",
    "alkene",
]


def find_smiles_column(df: pd.DataFrame) -> str:
    for col in df.columns:
        if col.lower() in _PREFERRED_SMILES_COLS and col.lower().endswith("smiles"):
            return col
    for col in df.columns:
        if col.lower() == "smiles":
            return col
    for col in df.columns:
        try:
            head = df[col].astype(str).head(50)
            ok_rate = head.apply(lambda s: Chem.MolFromSmiles(s) is not None).mean()
            if ok_rate > 0.8:
                return col
        except Exception:
            continue
    raise ValueError("Could not identify a SMILES column. Include a 'smiles' column or similar.")


def best_label_column(df: pd.DataFrame, smiles_col: str) -> Optional[str]:
    candidates = [c for c in df.columns if c != smiles_col]
    for p in ["substrate", "name", "alkene", "compound", "id", "label"]:
        for c in candidates:
            if c.lower() == p:
                return c
    for c in candidates:
        if df[c].dtype == object:
            return c
    return None


def find_table_csv(table_id: str, tables_dir: Path) -> Path:
    for cand in [
        tables_dir / f"{table_id}.csv",
        tables_dir / f"{table_id.replace('-', '_')}.csv",
        tables_dir / f"{table_id.replace('_', '-')}.csv",
    ]:
        if cand.exists():
            return cand
    raise FileNotFoundError(f"Could not locate CSV for table_id '{table_id}' in {tables_dir}")


# -----------------------------
# RDKit-safe sanitization and repairs
# -----------------------------

def sanitize_lenient(mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    """Sanitize while skipping kekulization so heteroaromatics do not fail.
    We keep aromatic flags intact so aromatic SMARTS like [c] continue to match.
    """
    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        return mol, None
    except Exception as e:
        return None, f"sanitize_failed: {e}"


def fix_overvalent_nitrogens(mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    """If any neutral nitrogens have valence 4, assign +1 formal charge and re-sanitize.
    This repairs products like N-acylpyridinium species that otherwise trip the
    "Explicit valence for atom N is greater than permitted" error.
    Returns (mol, err_msg).
    """
    try:
        rw = Chem.RWMol(mol)
        changed = False
        for a in rw.GetAtoms():
            if a.GetAtomicNum() != 7:
                continue
            if a.GetFormalCharge() == 0 and a.GetTotalValence() >= 4:
                a.SetFormalCharge(+1)
                changed = True
        if not changed:
            return mol, None
        fixed = rw.GetMol()
        fixed, err = sanitize_lenient(fixed)
        if fixed is None:
            return None, err or "sanitize_failed_after_fix"
        return fixed, None
    except Exception as e:
        return None, f"fix_overvalent_nitrogen_failed: {e}"


def parse_mol_from_smiles(smiles: str) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None, "MolFromSmiles returned None"
        return sanitize_lenient(mol)
    except Exception as e:
        return None, f"parse_failed: {e}"


def apply_template(template: str, reactant_mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    try:
        rxn = Reactions.ReactionFromSmarts(template)
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
        return prod, None
    except Exception as e:
        return None, f"reaction_failed: {e}"


def make_instance_reaction(reactant: Chem.Mol, product: Chem.Mol) -> Reactions.ChemicalReaction:
    rxn = Reactions.ChemicalReaction()
    rxn.AddReactantTemplate(reactant)
    rxn.AddProductTemplate(product)
    return rxn


def render_reaction_png(reactant: Chem.Mol, product: Chem.Mol, out_path: Path, sub_img_size: Tuple[int, int] = (250, 250)) -> None:
    # Ensure coordinates so depiction looks decent without kekulization
    for m in (reactant, product):
        try:
            rdDepictor.Compute2DCoords(m)
        except Exception:
            pass
    img = Draw.ReactionToImage(make_instance_reaction(reactant, product), subImgSize=sub_img_size)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    img.save(str(out_path))


# -----------------------------
# Main logic
# -----------------------------

def process_table(entry: Dict, tables_dir: Path, out_dir: Path, limit: int = 0, verbose: bool = True) -> None:
    table_id = entry["table_id"]
    template = entry["reaction_template"]

    csv_path = find_table_csv(table_id, tables_dir)
    if verbose:
        print(f"[INFO] Processing table '{table_id}' -> {csv_path}")
    df = pd.read_csv(csv_path)

    smiles_col = find_smiles_column(df)
    label_col = best_label_column(df, smiles_col)

    results: List[Dict] = []
    errors: List[Dict] = []

    for idx, row in df.iterrows():
        if limit and idx >= limit:
            break
        smiles_raw = str(row[smiles_col]).strip()
        if not smiles_raw or smiles_raw.lower() in {"nan", "none"}:
            continue

        label = str(row[label_col]) if label_col else smiles_raw
        slug = slugify(label)

        mol, err = parse_mol_from_smiles(smiles_raw)
        if mol is None:
            msg = err or "unknown_parse_error"
            errors.append({"row": idx, "label": label, "substrate_smiles": smiles_raw, "stage": "reactant_parse", "error": msg, "template": template})
            if verbose:
                print(f"[ERROR] {table_id} row {idx} ({label}) reactant_parse: {msg}  SMILES={smiles_raw}")
            continue

        prod, perr = apply_template(template, mol)
        if prod is None:
            msg = perr or "unknown_reaction_error"
            errors.append({"row": idx, "label": label, "substrate_smiles": smiles_raw, "stage": "reaction", "error": msg, "template": template})
            if verbose and msg != "no_products":
                print(f"[WARN] {table_id} row {idx} ({label}) reaction: {msg}")
            continue

        # Render
        png_path = out_dir / table_id / f"{idx:03d}_{slug}.png"
        try:
            render_reaction_png(mol, prod, png_path)
        except Exception:
            try:
                render_reaction_png(mol, prod, png_path, sub_img_size=(200, 200))
            except Exception as e2:
                errors.append({"row": idx, "label": label, "substrate_smiles": smiles_raw, "stage": "render", "error": f"{e2}", "template": template})
                if verbose:
                    print(f"[ERROR] {table_id} row {idx} ({label}) render: {e2}")
                continue

        results.append(
            {
                "row": idx,
                "label": label,
                "substrate_smiles": Chem.MolToSmiles(mol, kekuleSmiles=False),
                "product_smiles": Chem.MolToSmiles(prod, kekuleSmiles=False),
                "image": str(png_path.relative_to(out_dir)),
            }
        )

    out_tbl_dir = out_dir / table_id
    out_tbl_dir.mkdir(parents=True, exist_ok=True)
    if results:
        pd.DataFrame(results).to_csv(out_tbl_dir / "products.csv", index=False)
    if errors:
        err_df = pd.DataFrame(errors, columns=[
            "row", "label", "substrate_smiles", "stage", "error", "template"
        ])
        err_df.to_csv(out_tbl_dir / "errors.csv", index=False)
        if verbose:
            print(f"[INFO] {table_id}: wrote {len(errors)} error(s) to {out_tbl_dir / 'errors.csv'}")


def main(argv: Optional[Iterable[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Project reaction templates onto substrates and render reaction images.")
    ap.add_argument("--manifest", type=Path, required=True, help="Path to YAML manifest listing tables and reaction templates.")
    ap.add_argument("--tables-dir", type=Path, default=Path("."), help="Directory containing the CSVs (named <table_id>.csv).")
    ap.add_argument("--out-dir", type=Path, default=Path("analysis/reactivity_series/render"), help="Directory to write images and per-table products.csv")
    ap.add_argument("--limit", type=int, default=0, help="If >0, process only the first N rows of each table.")
    ap.add_argument("--quiet", action="store_true", help="Suppress info-level logs; still prints errors.")
    args = ap.parse_args(list(argv) if argv is not None else None)

    manifest = load_manifest(args.manifest)

    verbose = not args.quiet
    for entry in manifest:
        try:
            process_table(entry, args.tables_dir, args.out_dir, args.limit, verbose=verbose)
        except FileNotFoundError as e:
            print(f"[WARN] {e}")
        except Exception as e:
            print(f"[WARN] Failed to process {entry.get('table_id', '<unknown>')}: {e}")


if __name__ == "__main__":
    main()