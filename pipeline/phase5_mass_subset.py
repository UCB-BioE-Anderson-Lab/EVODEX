import csv
import os
import time
import sys
import json
import ast

from evodex.formula import calculate_exact_mass
from pipeline.version import __version__

csv.field_size_limit(sys.maxsize)

# Project root (.. from pipeline/)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
ERRORS_DIR = os.path.join(DATA_DIR, "errors")
EVODEX_DATA_DIR = os.path.join(BASE_DIR, "evodex", "data")

# Inputs (from Phase 3d)
EVODEX_F_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_f_phase3d_final.csv")
EVODEX_P_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_p_phase3d_final.csv")

# Outputs (pipeline)
EVODEX_M = os.path.join(PROCESSED_DIR, "evodex_m_phase5.csv")
EVODEX_M_SUBSET = os.path.join(PROCESSED_DIR, "evodex_m_mass_spec_subset_phase5.csv")


# Phase 5: Mass Spec Subset for EVODEX.2
# This phase generates EVODEX-M (full) and EVODEX-M_mass_spec_subset (filtered) tables.
#
# EVODEX-M contains exact mass differences for all formula differences in EVODEX-F.
# EVODEX-M_mass_spec_subset filters these to EVODEX-P reaction patterns compatible with
# single-species mass spec observations (at least one side of the reaction is a single fragment).
#
# EVODEX-M IDs are assigned as EVODEX.<__version__>-M1, M2, ...
# The published outputs are written into evodex/data:
#   - evodex/data/EVODEX-M.csv
#   - evodex/data/EVODEX-M_mass_spec_subset.csv


def ensure_directories():
    for path in [PROCESSED_DIR, ERRORS_DIR, EVODEX_DATA_DIR]:
        os.makedirs(path, exist_ok=True)


def print_stats(label, value):
    print(f"{label}: {value:,}")


def _parse_formula_diff_dictlike(formula_str: str) -> dict:
    """
    EVODEX-F formula_diff is stored in a dict-like form (stringified dict).
    Accept JSON dict strings or Python-literal dict strings.

    Examples:
      {"C": 1, "H": -2}
      {'C': 1, 'H': -2}
    """
    s = (formula_str or "").strip()
    if not s:
        return {}

    if not (s.startswith("{") and s.endswith("}")):
        raise ValueError(f"formula_diff is not dict-like: {s!r}")

    try:
        d = json.loads(s)
    except Exception:
        d = ast.literal_eval(s)

    if not isinstance(d, dict):
        raise ValueError(f"formula_diff did not parse to dict: {type(d)}")

    out = {}
    for k, v in d.items():
        if v is None:
            continue
        out[str(k)] = int(v)
    return out


def _mass_spec_compatible_p_smirks(smirks: str) -> bool:
    """
    Mass-spec compatibility filter:
      keep reactions where at least one side is a single species (no '.' separators).
    """
    left, right = smirks.split(">>")
    left_n = len([x for x in left.split(".") if x])
    right_n = len([x for x in right.split(".") if x])
    return (left_n == 1) or (right_n == 1)


def main():
    start_time = time.time()
    print("Phase 5 mass subset generation started...")

    ensure_directories()

    # Basic existence checks so a wrong path fails fast and loudly.
    if not os.path.exists(EVODEX_F_PHASE3D_FINAL):
        raise FileNotFoundError(f"File not found: {EVODEX_F_PHASE3D_FINAL}")
    if not os.path.exists(EVODEX_P_PHASE3D_FINAL):
        raise FileNotFoundError(f"File not found: {EVODEX_P_PHASE3D_FINAL}")

    total_f_rows = 0
    total_m_rows = 0
    total_subset_rows = 0
    errors_f = 0
    errors_p = 0

    first_f_errors = []
    first_p_errors = []

    print("Starting Phase 5: Mass spec subset generation...")

    # ------------------------------------------------------------------
    # Step 1: Generate EVODEX-M (full)
    # ------------------------------------------------------------------
    evodex_m_map = {}

    with open(EVODEX_F_PHASE3D_FINAL, "r", newline="") as infile:
        reader = csv.DictReader(infile)

        required = {"id", "sources", "formula_diff"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"EVODEX-F input missing columns: {sorted(missing)}; found: {reader.fieldnames}"
            )

        for row in reader:
            total_f_rows += 1
            try:
                orig_f_id = (row.get("id") or "").strip()
                if not orig_f_id:
                    raise ValueError("Missing EVODEX-F id")

                atom_diff = _parse_formula_diff_dictlike(row["formula_diff"])
                mass_diff = calculate_exact_mass(atom_diff)

                sources = set((row.get("sources") or "").split(","))
                sources = {s for s in sources if s}

                evodex_m_map[orig_f_id] = {
                    "mass": float(mass_diff),
                    "sources": sources,
                    "formula_diff": (row.get("formula_diff") or "").strip(),
                }
                total_m_rows += 1
            except Exception as e:
                errors_f += 1
                if len(first_f_errors) < 10:
                    first_f_errors.append(
                        (
                            (row.get("id") or "").strip(),
                            ((row.get("formula_diff") or "")[:200]),
                            str(e),
                        )
                    )

    if total_m_rows == 0:
        print("\nERROR: Parsed 0 EVODEX-M rows from EVODEX-F.")
        print_stats("Total EVODEX_F rows processed", total_f_rows)
        print_stats("Errors in EVODEX_F", errors_f)
        if first_f_errors:
            print("\nFirst EVODEX-F parse errors:")
            for fid, fstr, emsg in first_f_errors:
                print(f"  id={fid} formula_diff={fstr!r} error={emsg}")
        sys.exit(1)

    # Assign EVODEX-M IDs (rank by support count)
    sorted_m_entries = sorted(
        evodex_m_map.items(),
        key=lambda x: len(x[1]["sources"]),
        reverse=True,
    )

    m_id_map = {}
    for i, (orig_f_id, _) in enumerate(sorted_m_entries, start=1):
        m_id_map[orig_f_id] = f"EVODEX.{__version__}-M{i}"

    with open(EVODEX_M, "w", newline="") as outfile:
        writer = csv.DictWriter(
            outfile, fieldnames=["id", "mass", "sources", "formula_diff"]
        )
        writer.writeheader()
        for orig_f_id, data in sorted_m_entries:
            writer.writerow(
                {
                    "id": m_id_map[orig_f_id],
                    "mass": data["mass"],
                    "sources": ",".join(sorted(data["sources"])),
                    "formula_diff": data["formula_diff"],
                }
            )

    # ------------------------------------------------------------------
    # Step 2: Find mass-spec compatible EVODEX-P reactions
    # ------------------------------------------------------------------
    valid_p_ids = set()

    with open(EVODEX_P_PHASE3D_FINAL, "r", newline="") as infile:
        reader = csv.DictReader(infile)

        required = {"id", "smirks"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"EVODEX-P input missing columns: {sorted(missing)}; found: {reader.fieldnames}"
            )

        for row in reader:
            try:
                p_id = (row.get("id") or "").strip()
                sm = (row.get("smirks") or "").strip()
                if not p_id or not sm:
                    raise ValueError("Missing EVODEX-P id or smirks")

                if _mass_spec_compatible_p_smirks(sm):
                    valid_p_ids.add(p_id)
            except Exception as e:
                errors_p += 1
                if len(first_p_errors) < 10:
                    first_p_errors.append(
                        (
                            (row.get("id") or "").strip(),
                            ((row.get("smirks") or "")[:200]),
                            str(e),
                        )
                    )

    # ------------------------------------------------------------------
    # Step 3: Generate EVODEX-M_mass_spec_subset (filtered)
    # ------------------------------------------------------------------
    with open(EVODEX_M_SUBSET, "w", newline="") as outfile:
        writer = csv.DictWriter(
            outfile, fieldnames=["id", "mass", "sources", "formula_diff"]
        )
        writer.writeheader()

        for orig_f_id, data in evodex_m_map.items():
            intersecting_sources = data["sources"].intersection(valid_p_ids)
            if intersecting_sources:
                writer.writerow(
                    {
                        "id": m_id_map[orig_f_id],
                        "mass": data["mass"],
                        "sources": ",".join(sorted(intersecting_sources)),
                        "formula_diff": data["formula_diff"],
                    }
                )
                total_subset_rows += 1

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------
    print("\n=== Phase 5 Statistics ===")
    print_stats("Total EVODEX_F rows processed", total_f_rows)
    print_stats("Total EVODEX_M rows written", total_m_rows)
    print_stats("EVODEX_P valid (mass-spec compatible)", len(valid_p_ids))
    print_stats("Total EVODEX_M_mass_spec_subset rows written", total_subset_rows)
    print_stats("Errors in EVODEX_F", errors_f)
    print_stats("Errors in EVODEX_P", errors_p)

    if first_f_errors:
        print("\nFirst EVODEX-F parse errors:")
        for fid, fstr, emsg in first_f_errors:
            print(f"  id={fid} formula_diff={fstr!r} error={emsg}")

    if first_p_errors:
        print("\nFirst EVODEX-P parse errors:")
        for pid, sm, emsg in first_p_errors:
            print(f"  id={pid} smirks={sm!r} error={emsg}")

    print("Phase 5 complete: Mass tables written.")

    # ------------------------------------------------------------------
    # Publish to evodex/data
    # ------------------------------------------------------------------
    print("\nPublishing to evodex/data...")

    dst_m = os.path.join(EVODEX_DATA_DIR, "EVODEX-M.csv")
    with open(EVODEX_M, "r", newline="") as src_file, open(dst_m, "w", newline="") as dst_file:
        dst_file.write(src_file.read())
    print(f"Published EVODEX-M to {dst_m}")

    dst_m_subset = os.path.join(EVODEX_DATA_DIR, "EVODEX-M_mass_spec_subset.csv")
    with open(EVODEX_M_SUBSET, "r", newline="") as src_file, open(dst_m_subset, "w", newline="") as dst_file:
        dst_file.write(src_file.read())
    print(f"Published EVODEX-M_mass_spec_subset to {dst_m_subset}")

    elapsed_time = time.time() - start_time
    print(f"Phase 5 mass subset generation completed in {elapsed_time:.2f} seconds.")


if __name__ == "__main__":
    main()