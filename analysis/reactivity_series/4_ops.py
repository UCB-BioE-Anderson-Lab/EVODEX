#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Dict, Any

from evodex.operators import extract_operator_by_abstraction
from evodex.utils import reaction_hash

# Optional RdChiral dependency for radius=1, special_groups=True operators
try:
    from rdchiral import template_extractor as _rdchiral_te
except ImportError:
    _rdchiral_te = None

BASE = Path(__file__).resolve().parent
IN_DIR = BASE / "out/3_hmap"
OUT_DIR = BASE / "out/4_ops"

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
    """Attach EVODEX operators (A–E) to a single Stage-3 row.

    Expects a successful Stage-3 record that already contains:
      - "evodex_p_smirks": fully H-mapped EVODEX-P reaction SMIRKS

    On success, returns the row augmented with, for each abstraction level X in {A,B,C,D,E}:
      - "evodex_x_smirks": operator SMIRKS at abstraction level X
      - "evodex_x_hash": integer hash of that operator SMIRKS

    """

    smirks = row.get("evodex_p_smirks")
    if not smirks:
        raise ValueError("Missing evodex_p_smirks in Stage-3 row")

    out = dict(row)

    # Compute operators for all abstraction levels A–E
    for level in ("A", "B", "C", "D", "E"):
        op_smirks = extract_operator_by_abstraction(smirks, level)
        key = level.lower()
        out[f"evodex_{key}_smirks"] = op_smirks
        out[f"evodex_{key}_hash"] = reaction_hash(op_smirks)

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
        out["rdchiral_hash"] = reaction_hash(rdchiral_smirks)

    return out


def process_file(in_path: Path, out_path: Path) -> None:
    """Read a Stage-3 JSONL file and write the corresponding Stage-4 JSONL file.

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
        for line in fin:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)

            # If a previous stage marked this row as an error, just propagate it
            if "error" in rec:
                fout.write(json.dumps(rec, ensure_ascii=False, indent=2) + "\n")
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
                fout.write(json.dumps(err_rec, ensure_ascii=False, indent=2) + "\n")
            else:
                fout.write(json.dumps(out_rec, ensure_ascii=False, indent=2) + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    table_filter = set(ONLY_TABLES)

    input_files = sorted(p for p in IN_DIR.glob("*.jsonl"))
    if table_filter:
        input_files = [p for p in input_files if p.stem in table_filter]

    if not input_files:
        print(f"[4_ops] No input files found in {IN_DIR} (filter={sorted(table_filter)})")
        return

    print(f"[4_ops] Writing Stage-4 operator files to {OUT_DIR}")

    for in_path in input_files:
        out_path = OUT_DIR / in_path.name
        print(f"[4_ops] Processing {in_path.name} -> {out_path.name}")
        process_file(in_path, out_path)


if __name__ == "__main__":
    main()
