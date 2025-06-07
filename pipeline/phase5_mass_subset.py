import csv
import os
from collections import defaultdict
from evodex.formula import calculate_exact_mass
from pipeline.config import load_paths
from pipeline.version import __version__

# Phase 5: Mass Spec Subset
# This phase generates EVODEX_M and EVODEX_M_SUBSET tables.
# EVODEX_M contains exact mass differences for all formula differences in EVODEX_F.
# EVODEX_M_SUBSET filters these to reaction patterns (EVODEX_P) compatible with single-species
# mass spec observations (i.e., at least one side of the reaction is a single fragment).

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def print_stats(label, value):
    print(f"{label}: {value:,}")

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    total_f_rows = 0
    total_m_rows = 0
    total_subset_rows = 0
    errors_f = 0
    errors_p = 0

    print("Starting Phase 5: Mass spec subset generation...")
    
    # Step 1: generate evodex_m
    evodex_m_map = {}
    with open(paths['evodex_f_phase3c_final'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            total_f_rows += 1
            try:
                formula_diff = eval(row['formula'])
                mass_diff = calculate_exact_mass(formula_diff)
                evodex_m_map[row['id']] = {
                    'mass': mass_diff,
                    'sources': set(row['sources'].split(',')),
                    'formula': row['formula']
                }
                total_m_rows += 1
            except Exception:
                errors_f += 1
                continue

    # Step 1b: assign EVODEX_M IDs
    sorted_m_entries = sorted(evodex_m_map.items(), key=lambda x: len(x[1]['sources']), reverse=True)
    m_id_map = {}
    for i, (orig_id, data) in enumerate(sorted_m_entries, start=1):
        assigned_id = f"EVODEX.1-M{i}"
        m_id_map[orig_id] = assigned_id

    # Write EVODEX_M with new IDs
    with open(paths['evodex_m'], 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'mass', 'sources', 'formula'])
        writer.writeheader()
        for orig_id, data in sorted_m_entries:
            assigned_id = m_id_map[orig_id]
            writer.writerow({
                'id': assigned_id,
                'mass': data['mass'],
                'sources': ','.join(sorted(data['sources'])),
                'formula': data['formula']
            })

    # Step 2: filter evodex_p for mass-spec compatible ones
    valid_p_ids = set()
    with open(paths['evodex_p_phase3c_final'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            try:
                left, right = row['smirks'].split('>>')
                if len(left.split('.')) == 1 or len(right.split('.')) == 1:
                    valid_p_ids.add(row['id'])
            except Exception:
                errors_p += 1
                continue

    # Step 3: generate evodex_m_subset
    with open(paths['evodex_m_subset'], 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'mass', 'sources'])
        writer.writeheader()
        for orig_id, data in evodex_m_map.items():
            intersecting_sources = data['sources'].intersection(valid_p_ids)
            if intersecting_sources:
                assigned_id = m_id_map[orig_id]
                writer.writerow({
                    'id': assigned_id,
                    'mass': data['mass'],
                    'sources': ','.join(sorted(intersecting_sources))
                })
                total_subset_rows += 1

    print("\n=== Phase 5 Statistics ===")
    print_stats("Total EVODEX_F rows processed", total_f_rows)
    print_stats("Total EVODEX_M rows written", total_m_rows)
    print_stats("EVODEX_P valid (mass-spec compatible)", len(valid_p_ids))
    print_stats("Total EVODEX_M_SUBSET rows written", total_subset_rows)
    print_stats("Errors in EVODEX_F", errors_f)
    print_stats("Errors in EVODEX_P", errors_p)

    print("Phase 5 complete: Mass subset written.")

    # Publish to evodex/data
    print("\nPublishing to evodex/data...")

    # EVODEX-M
    dst_m = os.path.join('evodex', 'data', 'EVODEX-M_unique_masses.csv')
    with open(paths['evodex_m'], 'r') as src_file, open(dst_m, 'w', newline='') as dst_file:
        dst_file.write(src_file.read())
    print(f"Published EVODEX-M to {dst_m}")

    # EVODEX-M_SUBSET
    dst_m_subset = os.path.join('evodex', 'data', 'EVODEX-M_mass_spec_subset.csv')
    with open(paths['evodex_m_subset'], 'r') as src_file, open(dst_m_subset, 'w', newline='') as dst_file:
        dst_file.write(src_file.read())
    print(f"Published EVODEX-M_SUBSET to {dst_m_subset}")

if __name__ == "__main__":
    main()