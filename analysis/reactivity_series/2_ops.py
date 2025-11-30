#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Dict, Any

from evodex.operators import extract_operator_by_abstraction

# Optional RdChiral dependency for radius=1, special_groups=True operators
try:
    from rdchiral import template_extractor as _rdchiral_te
except ImportError:
    _rdchiral_te = None

BASE = Path(__file__).resolve().parent
IN_DIR = BASE / "out/1_project"
OUT_DIR = BASE / "out/2_ops"
LOG_DIR = BASE / "out/_logs"

# If nonempty, only process tables whose file stem (for example "alcohol_pka") is listed here
ONLY_TABLES: Iterable[str] = []

EVODEX_ABSTRACTION_LEVELS = ("A", "B", "C", "D", "E")


def _rdchiral_operator_from_smirks(forward_smirks: str) -> str | None:
    """Return a RdChiral operator SMIRKS for a forward reaction SMIRKS.

    The operator uses radius=1 for reactants, radius=0 for products and
    special_groups=True. RdChiral templates are retrosynthetic, so we pass
    the forward reaction to RdChiral with reactants and products swapped
    and force the desired radii via get_fragments_for_changed_atoms.
    """
    if _rdchiral_te is None:
        return None

    if ">>" not in forward_smirks:
        raise ValueError(f"SMIRKS missing '>>' separator: {forward_smirks}")

    lhs, rhs = forward_smirks.split(">>", 1)
    rd_reaction = {
        "reactants": rhs,
        "products": lhs,
        "_id": "evodex",
    }

    te = _rdchiral_te

    # Snapshot rdchiral global configuration
    orig_VERBOSE = getattr(te, "VERBOSE", False)
    orig_USE_STEREOCHEMISTRY = getattr(te, "USE_STEREOCHEMISTRY", True)
    orig_MAX_UNMAPPED = getattr(te, "MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS", 100)
    orig_INCLUDE_ALL_UNMAPPED = getattr(te, "INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS", True)
    orig_get_fragments = te.get_fragments_for_changed_atoms
    orig_get_special_groups = te.get_special_groups

    # Configure radius and other settings
    te.VERBOSE = False
    te.USE_STEREOCHEMISTRY = True
    te.MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS = 100
    te.INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS = True

    # Wrap get_fragments_for_changed_atoms to enforce radii
    def _gffca_forced(mols, changed_atom_tags, radius=0, category="reactants", expansion=None):
        if expansion is None:
            expansion = []
        forced_radius = 1 if category == "reactants" else 0
        return orig_get_fragments(
            mols,
            changed_atom_tags,
            radius=forced_radius,
            category=category,
            expansion=expansion,
        )

    te.get_fragments_for_changed_atoms = _gffca_forced

    try:
        tpl = te.extract_from_reaction(rd_reaction)
    finally:
        # Restore original globals to avoid side effects elsewhere
        te.VERBOSE = orig_VERBOSE
        te.USE_STEREOCHEMISTRY = orig_USE_STEREOCHEMISTRY
        te.MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS = orig_MAX_UNMAPPED
        te.INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS = orig_INCLUDE_ALL_UNMAPPED
        te.get_fragments_for_changed_atoms = orig_get_fragments
        te.get_special_groups = orig_get_special_groups

    # RdChiral returns None or a dict without reaction_smarts on failure
    if not tpl or "reaction_smarts" not in tpl:
        return None
    return tpl["reaction_smarts"]


def process_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Attach EVODEX operators (A to E) and optional RdChiral operator to one row.

    Expects a successful record from 1_project that already contains
    "mapped_substrate_smiles" and "mapped_product_smiles". These are combined
    into a forward EVODEX-P SMIRKS stored as "evodex_p_smirks".

    On success, returns a copy of the row with the following extra fields:
      - "evodex_p_smirks"
      - "evodex_a_smirks", ..., "evodex_e_smirks"
      - "rdchiral_smirks" (possibly None)
      - "rdchiral_error" only if RdChiral was unavailable or failed
    """
    mapped_sub = row.get("mapped_substrate_smiles")
    mapped_prod = row.get("mapped_product_smiles")

    if not mapped_sub or not mapped_prod:
        raise ValueError("Missing mapped_substrate_smiles or mapped_product_smiles in projection row")

    smirks = f"{mapped_sub}>>{mapped_prod}"

    out = dict(row)
    out["evodeX_p_smirks"] = smirks  # keep spelling consistent with prior stages if applicable
    out["evodex_p_smirks"] = smirks

    for level in EVODEX_ABSTRACTION_LEVELS:
        op_smirks = extract_operator_by_abstraction(smirks, level)
        key = level.lower()
        out[f"evodex_{key}_smirks"] = op_smirks

    rdchiral_smirks = None
    if _rdchiral_te is None:
        out["rdchiral_error"] = (
            "rdchiral not installed; install with "
            "'pip install git+https://github.com/connorcoley/rdchiral.git'"
        )
    else:
        try:
            rdchiral_smirks = _rdchiral_operator_from_smirks(smirks)
            if not rdchiral_smirks:
                out["rdchiral_error"] = "rdchiral template_extractor returned no template"
        except Exception as e:
            rdchiral_smirks = None
            out["rdchiral_error"] = f"{type(e).__name__}: {e}"

    out["rdchiral_smirks"] = rdchiral_smirks

    return out


def _make_error_record(stage: str, error: str, rec: Dict[str, Any]) -> Dict[str, Any]:
    """Build a compact error record compatible with earlier pipeline stages."""
    err_rec: Dict[str, Any] = {
        "stage": stage,
        "error": error,
    }
    for key in ("table_id", "row", "label"):
        if key in rec:
            err_rec[key] = rec[key]
    return err_rec


def _log_rdchiral_error(out_rec: Dict[str, Any], err_f) -> None:
    """Write a nonfatal RdChiral error record if present in the row."""
    if "rdchiral_error" not in out_rec:
        return
    rd_err = _make_error_record("operator_extraction_rdchiral", out_rec["rdchiral_error"], out_rec)
    err_f.write(json.dumps(rd_err, ensure_ascii=False) + "\n")


def process_file(in_path: Path, out_path: Path, err_f) -> None:
    """Convert a Stage 1 projection JSONL file to a Stage 2 operator JSONL file.

    Error handling is row-local: if operator extraction fails for a row we emit
    a compact error record with stage set to "operator_extraction". Successful
    rows are written with EVODEX operator fields added.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with in_path.open("r", encoding="utf-8") as fin, out_path.open("w", encoding="utf-8") as fout:
        buffer: list[str] = []

        def _flush_buffer():
            nonlocal buffer
            if not buffer:
                return None
            text = "\n".join(buffer)
            buffer = []
            return json.loads(text)

        for raw_line in fin:
            line = raw_line.rstrip("\n")
            if not line.strip():
                continue

            stripped = line.strip()
            if not buffer and stripped.startswith("{") and stripped.endswith("}"):
                rec = json.loads(stripped)
            else:
                buffer.append(line)
                if stripped.endswith("}"):
                    rec = _flush_buffer()
                else:
                    continue

            if rec is None:
                continue

            if "error" in rec:
                fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
                continue

            try:
                out_rec = process_row(rec)
                _log_rdchiral_error(out_rec, err_f)
            except Exception as e:
                err_rec = _make_error_record("operator_extraction", str(e), rec)
                err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                fout.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
            else:
                fout.write(json.dumps(out_rec, ensure_ascii=False) + "\n")

        if buffer:
            rec = _flush_buffer()
            if rec is not None:
                if "error" in rec:
                    fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
                else:
                    try:
                        out_rec = process_row(rec)
                    except Exception as e:
                        err_rec = _make_error_record("operator_extraction", str(e), rec)
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        fout.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                    else:
                        _log_rdchiral_error(out_rec, err_f)
                        fout.write(json.dumps(out_rec, ensure_ascii=False) + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    errors_path = LOG_DIR / "2_ops.errors.jsonl"

    table_filter = set(ONLY_TABLES)

    input_files = sorted(p for p in IN_DIR.glob("*.jsonl"))
    if table_filter:
        input_files = [p for p in input_files if p.stem in table_filter]

    if not input_files:
        print(f"[2_ops] No input files found in {IN_DIR} (filter={sorted(table_filter)})")
        return

    print(f"[2_ops] Writing Stage 4 operator files to {OUT_DIR}")

    with open(errors_path, "w", encoding="utf-8") as err_f:
        for in_path in input_files:
            out_path = OUT_DIR / in_path.name
            print(f"[2_ops] Processing {in_path.name} -> {out_path.name}")
            process_file(in_path, out_path, err_f)


if __name__ == "__main__":
    main()
