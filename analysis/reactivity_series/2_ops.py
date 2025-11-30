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

# If nonempty, only process tables whose file stem (e.g. "alcohol_pka") is listed here
ONLY_TABLES: Iterable[str] = []

# Helper to compute RdChiral operator from a forward SMIRKS string
def _rdchiral_operator_from_smirks(forward_smirks: str) -> str | None:
    """Return RdChiral radius=1 (reactants), 0 (products), special_groups=True operator
    SMIRKS for a forward reaction SMIRKS.

    RdChiral templates are retrosynthetic by design, so we pass the forward
    reaction to RdChiral with reactants/products swapped (RHS as reactants,
    LHS as products) and force the desired radius via the internal
    `get_fragments_for_changed_atoms` helper.
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

    # Snapshot rdchiral's global configuration
    orig_VERBOSE = getattr(te, "VERBOSE", False)
    orig_USE_STEREOCHEMISTRY = getattr(te, "USE_STEREOCHEMISTRY", True)
    orig_MAX_UNMAPPED = getattr(te, "MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS", 5)
    orig_INCLUDE_ALL_UNMAPPED = getattr(te, "INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS", True)
    orig_get_fragments = te.get_fragments_for_changed_atoms
    orig_get_special_groups = te.get_special_groups

    # Configure for our abstraction: radius=1 (reactants), 0 (products),
    # special_groups=True, stereo on, and a modest unmapped-atom limit.
    te.VERBOSE = False
    te.USE_STEREOCHEMISTRY = True
    te.MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS = 5
    te.INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS = True  # or False if we want a tighter core

    # Keep special groups enabled (default behavior) by leaving get_special_groups as-is.

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
    """Attach EVODEX operators (A–E) to a single Stage-1 projection row.

    Expects a successful record from 1_project that already contains:
      - "mapped_substrate_smiles": atom-mapped substrate SMILES
      - "mapped_product_smiles": atom-mapped product SMILES

    These are combined into a forward EVODEX-P SMIRKS of the form
    "evodex_p_smirks" = "{mapped_substrate_smiles}>>{mapped_product_smiles}".

    On success, returns the row augmented with:
      - "evodex_p_smirks": forward EVODEX-P reaction SMIRKS
      and for each abstraction level X in {A,B,C,D,E}:
      - "evodex_x_smirks": operator SMIRKS at abstraction level X

    If RdChiral is available, the row may also be augmented with:
      - "rdchiral_smirks": RdChiral-style operator SMIRKS
    """
    mapped_sub = row.get("mapped_substrate_smiles")
    mapped_prod = row.get("mapped_product_smiles")

    if not mapped_sub or not mapped_prod:
        raise ValueError("Missing mapped_substrate_smiles or mapped_product_smiles in projection row")

    # Forward EVODEX-P SMIRKS constructed directly from projection output
    smirks = f"{mapped_sub}>>{mapped_prod}"

    out = dict(row)
    out["evodex_p_smirks"] = smirks

    # Compute operators for all abstraction levels A–E
    for level in ("A", "B", "C", "D", "E"):
        op_smirks = extract_operator_by_abstraction(smirks, level)
        key = level.lower()
        out[f"evodex_{key}_smirks"] = op_smirks

    # Optionally compute a RdChiral operator (radius=1, special_groups=True)
    # from the same fully H-mapped EVODEX-P SMIRKS.
    rdchiral_smirks = None
    if _rdchiral_te is None:
        # Soft failure: record that RdChiral is unavailable but do not abort the row
        out["rdchiral_error"] = (
            "rdchiral not installed; install with "
            "'pip install git+https://github.com/connorcoley/rdchiral.git'"
        )
    else:
        try:
            rdchiral_smirks = _rdchiral_operator_from_smirks(smirks)
        except Exception as e:
            # Keep the EVODEX fields even if RdChiral fails
            out["rdchiral_error"] = f"{type(e).__name__}: {e}"

    if rdchiral_smirks:
        out["rdchiral_smirks"] = rdchiral_smirks

    return out


def process_file(in_path: Path, out_path: Path, err_f) -> None:
    """Read a Stage-1 projection JSONL file and write the corresponding Stage-2 operator JSONL file.

    Error handling is row-local: if operator extraction fails for a row
    we emit a compact error record of the form used by earlier stages:

      {
        "table_id": ..., "row": ..., "label": ...,
        "stage": "operator_extraction",
        "error": "...stringified exception..."
      }

    Successful rows are passed through with new EVODEX-E fields added.
    """

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with in_path.open("r") as fin, out_path.open("w") as fout:
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
                # skip blank lines
                continue

            # Fast path: a single-line JSON object
            stripped = line.strip()
            if not buffer and stripped.startswith("{") and stripped.endswith("}"):
                rec = json.loads(stripped)
            else:
                # Accumulate lines until we hit the end of an object (line ending with "}")
                buffer.append(line)
                if stripped.endswith("}"):
                    rec = _flush_buffer()
                else:
                    # Need more lines to complete this object
                    continue

            if rec is None:
                # Should not happen, but guard anyway
                continue

            # If a previous stage marked this row as an error, just propagate it
            if "error" in rec:
                fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
                continue

            try:
                out_rec = process_row(rec)
            except Exception as e:
                # Fall back to the minimal error schema we have from earlier stages
                err_rec: Dict[str, Any] = {
                    "stage": "operator_extraction",
                    "error": str(e),
                }
                # Preserve basic identifiers if present
                for key in ("table_id", "row", "label"):
                    if key in rec:
                        err_rec[key] = rec[key]
                err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                fout.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
            else:
                fout.write(json.dumps(out_rec, ensure_ascii=False) + "\n")

        # In case there is a final buffered object without trailing newline
        if buffer:
            rec = _flush_buffer()
            if rec is not None:
                if "error" in rec:
                    fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
                else:
                    try:
                        out_rec = process_row(rec)
                    except Exception as e:
                        err_rec = {
                            "stage": "operator_extraction",
                            "error": str(e),
                        }
                        for key in ("table_id", "row", "label"):
                            if key in rec:
                                err_rec[key] = rec[key]
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        fout.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                    else:
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

    print(f"[2_ops] Writing Stage-4 operator files to {OUT_DIR}")

    with open(errors_path, "w", encoding="utf-8") as err_f:
        for in_path in input_files:
            out_path = OUT_DIR / in_path.name
            print(f"[2_ops] Processing {in_path.name} -> {out_path.name}")
            process_file(in_path, out_path, err_f)


if __name__ == "__main__":
    main()
