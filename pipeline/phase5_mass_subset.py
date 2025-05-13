

import csv
import os
from collections import defaultdict
from evodex.formula import calculate_exact_mass
from pipeline.config import load_paths
from pipeline.version import __version__

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    print("Starting Phase 5: Mass spec subset generation...")
    
    # Step 1: generate evodex_m
    evodex_m_map = {}
    with open(paths['evodex_f'], 'r') as infile, open(paths['evodex_m'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=['id', 'mass', 'sources', 'formula'])
        writer.writeheader()
        for row in reader:
            try:
                formula_diff = eval(row['formula'])
                mass_diff = calculate_exact_mass(formula_diff)
                evodex_m_map[row['id']] = {
                    'mass': mass_diff,
                    'sources': set(row['sources'].split(',')),
                    'formula': row['formula']
                }
                writer.writerow({
                    'id': row['id'],
                    'mass': mass_diff,
                    'sources': row['sources'],
                    'formula': row['formula']
                })
            except Exception:
                continue

    # Step 2: filter evodex_p for mass-spec compatible ones
    valid_p_ids = set()
    with open(paths['evodex_p'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            try:
                left, right = row['smirks'].split('>>')
                if len(left.split('.')) == 1 or len(right.split('.')) == 1:
                    valid_p_ids.add(row['id'])
            except Exception:
                continue

    # Step 3: generate evodex_m_subset
    with open(paths['evodex_m_subset'], 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'mass', 'sources'])
        writer.writeheader()
        for m_id, m_data in evodex_m_map.items():
            intersecting_sources = m_data['sources'].intersection(valid_p_ids)
            if intersecting_sources:
                writer.writerow({
                    'id': m_id,
                    'mass': m_data['mass'],
                    'sources': ','.join(sorted(intersecting_sources))
                })

    print("Phase 5 complete: Mass subset written.")

if __name__ == "__main__":
    main()