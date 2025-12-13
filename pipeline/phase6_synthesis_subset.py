import csv
import os
import shutil
import time
from collections import defaultdict
from pipeline.version import __version__
import sys
csv.field_size_limit(sys.maxsize)

# Project root (.. from pipeline/)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
EVODEX_DATA_DIR = os.path.join(BASE_DIR, "evodex", "data")

# Inputs from Phase 3d / Phase 4
EVODEX_D_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_d_complete.csv")
EVODEX_P_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_p_phase3d_final.csv")

# Synthesis subset output (write directly into evodex/data)
EVODEX_D_SYNTHESIS = os.path.join(EVODEX_DATA_DIR, "EVODEX-D_synthesis_subset.csv")

# Phase 6: Synthesis Subset (EVODEX-D)
# This phase filters EVODEX-D operators to a subset usable for synthesis algorithms.
# It iterates through EVODEX-P reactions and filters based on single-species compatibility.
# The remaining EVODEX-P IDs are used to select associated EVODEX-D operators, which are then output as the synthesis subset.
#
# Note: The 'sources' column in the synthesis subset reflects only EVODEX-P entries
# that passed the compatibility filter. Therefore, the source counts and order of EVODEX-D entries
# may differ from the full EVODEX-D table produced in Phase 4 (which used all EVODEX-P sources).
# EVODEX-D IDs remain unchanged and consistent across phases.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)

def _is_single_species_side(side: str) -> bool:
    parts = [p for p in (side or "").split(".") if p]
    return len(parts) == 1

def _is_synthesis_compatible_p(p_smirks: str) -> bool:
    # Single-species compatibility filter based on EVODEX-P reaction SMIRKS.
    # We keep reactions where at least one side is a single species (single fragment),
    # which is compatible with common single-species observation / synthesis projections.
    left, right = p_smirks.split(">>")
    return _is_single_species_side(left) or _is_single_species_side(right)

def main():
    start_time = time.time()
    print("Phase 6 synthesis subset generation started...")
    print("Starting Phase 6: EVODEX-D synthesis subset generation (compatibility filtering)...")

    ensure_directories({"EVODEX_D_SYNTHESIS": EVODEX_D_SYNTHESIS})

    if not os.path.exists(EVODEX_P_PHASE3D_FINAL):
        raise FileNotFoundError(f"File not found: {EVODEX_P_PHASE3D_FINAL}")
    if not os.path.exists(EVODEX_D_PHASE3D_FINAL):
        raise FileNotFoundError(f"File not found: {EVODEX_D_PHASE3D_FINAL}")

    evodex_p_map = {}
    evodex_d_map = {}
    evodex_d_full_map = {}

    # Load EVODEX-P reactions (processed Phase 3d final)
    with open(EVODEX_P_PHASE3D_FINAL, "r", newline="") as p_file:
        p_reader = csv.DictReader(p_file)
        for row in p_reader:
            pid = (row.get("id") or "").strip()
            sm = (row.get("smirks") or "").strip()
            if pid and sm:
                evodex_p_map[pid] = sm

    # Load EVODEX-D operators mapping from EVODEX-P IDs (processed D-complete)
    with open(EVODEX_D_PHASE3D_FINAL, "r", newline="") as d_file:
        d_reader = csv.DictReader(d_file)
        for row in d_reader:
            did = (row.get("id") or "").strip()
            d_smirks = (row.get("smirks") or "").strip()
            sources = [s for s in (row.get("sources") or "").split(",") if s]

            if not did or not d_smirks:
                continue

            evodex_d_full_map[did] = {"smirks": d_smirks, "sources": sources}
            for p_id in sources:
                evodex_d_map.setdefault(p_id, []).append(did)

    evodex_d_subset = set()
    evodex_d_sources = {}
    p_parse_errors = 0

    # Process EVODEX-P reactions, filtering based on single-species compatibility
    for p_id, p_smirks in evodex_p_map.items():
        try:
            if not _is_synthesis_compatible_p(p_smirks):
                continue
        except Exception:
            p_parse_errors += 1
            continue

        # Add all D operators supported by this compatible P-reaction
        for did in evodex_d_map.get(p_id, []):
            d_smirks = evodex_d_full_map[did]["smirks"]

            # Keep the old operator-side guard to ensure the operator itself is single-species
            # (avoid multi-component operators for synthesis projection).
            if "." in d_smirks:
                continue

            evodex_d_subset.add(did)
            evodex_d_sources.setdefault(did, set()).add(p_id)

    # Write filtered EVODEX-D synthesis subset
    with open(EVODEX_D_SYNTHESIS, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["id", "smirks", "sources"])
        writer.writeheader()

        # Stable ordering: most supported first, then ID
        for did in sorted(
            evodex_d_subset,
            key=lambda x: (-len(evodex_d_sources.get(x, set())), x),
        ):
            writer.writerow(
                {
                    "id": did,
                    "smirks": evodex_d_full_map[did]["smirks"],
                    "sources": ",".join(sorted(evodex_d_sources.get(did, set()))),
                }
            )

    # Summary statistics
    total_evodex_p = len(evodex_p_map)
    total_evodex_p_filtered = len(set().union(*evodex_d_sources.values())) if evodex_d_sources else 0
    total_evodex_d_written = len(evodex_d_subset)

    print("Phase 6 Statistics:")
    print(f"  Total EVODEX-P entries processed: {total_evodex_p}")
    print(f"  Total EVODEX-P entries after compatibility filtering: {total_evodex_p_filtered}")
    print(f"  Total EVODEX-D synthesis operators written: {total_evodex_d_written}")
    if p_parse_errors:
        print(f"  EVODEX-P SMIRKS parse errors: {p_parse_errors}")

    print("Phase 6 complete: Synthesis subset written.")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 6 synthesis subset generation completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()