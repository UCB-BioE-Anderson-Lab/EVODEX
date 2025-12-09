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

# Inputs from Phase 3d
EVODEX_D_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_d_complete.csv")
EVODEX_P_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_p_phase3d_final.csv")

# Synthesis subset output (processed)
EVODEX_D_SYNTHESIS = os.path.join(PROCESSED_DIR, "EVODEX-D_synthesis_subset.csv")

# Phase 6: Synthesis Subset
# This phase filters EVODEX-E operators to a subset usable for synthesis algorithms.
# It iterates through EVODEX-P reactions and filters based on single-substrate/single-product compatibility.
# The remaining EVODEX-P IDs are used to find associated EVODEX-E operators, which are then output as the synthesis subset.
# Note: The 'sources' column in the synthesis subset reflects only EVODEX-P entries
# that passed the compatibility filter. Therefore, the source counts and order of EVODEX-E entries
# may differ from the full EVODEX-E table produced in Phase 3 (which used all EVODEX-P sources).
# EVODEX-E IDs remain unchanged and consistent across phases.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    start_time = time.time()
    print("Phase 6 synthesis subset generation started...")

    print("Starting Phase 6: Synthesis subset generation (compatibility filtering)...")

    evodex_p_map = {}
    evodex_d_map = {}
    evodex_d_full_map = {}

    # Load EVODEX-P reactions
    with open('evodex/data/EVODEX-P_partial_reactions.csv', 'r') as p_file:
        p_reader = csv.DictReader(p_file)
        for row in p_reader:
            evodex_p_map[row['id']] = row['smirks']

    # Load EVODEX-D operators mapping from EVODEX-P IDs
    with open(EVODEX_D_PHASE3D_FINAL, 'r') as e_file:
        e_reader = csv.DictReader(e_file)
        for row in e_reader:
            sources = row['sources'].split(',')
            evodex_d_full_map[row['id']] = {
                'smirks': row['smirks'],
                'sources': sources
            }
            for source in sources:
                evodex_d_map.setdefault(source, []).append(row['id'])

    evodex_d_subset = set()
    evodex_d_sources = {}

    # Process EVODEX-P reactions, filtering based on single-substrate/single-product compatibility
    for p_id, smirks in evodex_p_map.items():
        evodex_d_ids = evodex_d_map.get(p_id)
        if evodex_d_ids:
            for did in evodex_d_ids:
                d_smirks = evodex_d_full_map[did]['smirks']
                if '.' in d_smirks:
                    continue
                evodex_d_subset.add(did)
                evodex_d_sources.setdefault(did, set()).add(p_id)

    # Write filtered EVODEX-D synthesis subset
    with open(EVODEX_D_SYNTHESIS, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'smirks', 'sources'])
        writer.writeheader()
        for did in sorted(evodex_d_subset):
            row_data = evodex_d_full_map[did]
            writer.writerow({
                'id': did,
                'sources': ','.join(sorted(evodex_d_sources[did])),
                'smirks': row_data['smirks']
            })

    # Summary statistics
    total_evode_p = len(evodex_p_map)
    total_evode_p_filtered = len(set().union(*evodex_d_sources.values())) if evodex_d_sources else 0
    total_evode_d_written = len(evodex_d_subset)

    print(f"Phase 6 Statistics:")
    print(f"  Total EVODEX-P entries processed: {total_evode_p}")
    print(f"  Total EVODEX-P entries after compatibility filtering: {total_evode_p_filtered}")
    print(f"  Total EVODEX-D synthesis operators written: {total_evode_d_written}")

    print("Phase 6 complete: Synthesis subset written.")

    # === Publish to evodex/data/ ===
    dst_d_synthesis = os.path.join('evodex', 'data', 'EVODEX-D_synthesis_subset.csv')
    shutil.copyfile(EVODEX_D_SYNTHESIS, dst_d_synthesis)
    print(f"Published synthesis subset to {dst_d_synthesis}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 6 synthesis subset generation completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()