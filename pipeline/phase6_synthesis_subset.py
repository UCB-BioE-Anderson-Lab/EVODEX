import csv
import os
from collections import defaultdict
from evodex.decofactor import contains_cofactor
from pipeline.config import load_paths
from pipeline.version import __version__

# Phase 6: Synthesis Subset
# This phase filters EVODEX-E operators to a subset usable for synthesis algorithms.
# It iterates through EVODEX-P reactions, skipping any reactions
# that involve cofactors as substrates or products. The remaining EVODEX-P IDs
# are used to find associated EVODEX-E operators, which are then output as the synthesis subset.
# Note: The 'sources' column in the synthesis subset reflects only EVODEX-P entries
# that passed the cofactor filter. Therefore, the source counts and order of EVODEX-E entries
# may differ from the full EVODEX-E table produced in Phase 3 (which used all EVODEX-P sources).
# EVODEX-E IDs remain unchanged and consistent across phases.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    print("Starting Phase 6: Synthesis subset generation (cofactor filtering)...")

    evodex_p_map = {}
    evodex_e_map = {}
    evodex_e_full_map = {}

    # Load EVODEX-P reactions
    with open(paths['evodex_p'], 'r') as p_file:
        p_reader = csv.DictReader(p_file)
        for row in p_reader:
            evodex_p_map[row['id']] = row['smirks']

    # Load EVODEX-E operators mapping from EVODEX-P IDs
    with open(paths['evodex_e'], 'r') as e_file:
        e_reader = csv.DictReader(e_file)
        for row in e_reader:
            sources = row['sources'].split(',')
            evodex_e_full_map[row['id']] = {
                'smirks': row['smirks'],
                'sources': sources
            }
            for source in sources:
                evodex_e_map.setdefault(source, []).append(row['id'])

    evodex_e_subset = set()
    evodex_e_sources = {}

    # Process EVODEX-P reactions, skipping those with cofactors in substrates or products
    for p_id, smirks in evodex_p_map.items():
        try:
            # Skip if reaction contains any cofactor
            if contains_cofactor(smirks):
                continue

            # If passed filtering, find associated EVODEX-E IDs
            evodex_e_ids = evodex_e_map.get(p_id)
            if evodex_e_ids:
                for eid in evodex_e_ids:
                    evodex_e_subset.add(eid)
                    evodex_e_sources.setdefault(eid, set()).add(p_id)

        except Exception:
            # Skip malformed entries or parsing errors
            continue

    # Write filtered EVODEX-E synthesis subset
    with open(paths['evodex_e_synthesis'], 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'smirks', 'sources'])
        writer.writeheader()
        for eid in sorted(evodex_e_subset):
            row_data = evodex_e_full_map[eid]
            writer.writerow({
                'id': eid,
                'sources': ','.join(sorted(evodex_e_sources[eid])),
                'smirks': row_data['smirks']
            })

    # Summary statistics
    total_evode_p = len(evodex_p_map)
    total_evode_p_filtered = len(set().union(*evodex_e_sources.values())) if evodex_e_sources else 0
    total_evode_e_written = len(evodex_e_subset)

    print(f"Phase 6 Statistics:")
    print(f"  Total EVODEX-P entries processed: {total_evode_p}")
    print(f"  Total EVODEX-P entries after cofactor filtering: {total_evode_p_filtered}")
    print(f"  Total EVODEX-E synthesis operators written: {total_evode_e_written}")

    print("Phase 6 complete: Synthesis subset written.")

if __name__ == "__main__":
    main()