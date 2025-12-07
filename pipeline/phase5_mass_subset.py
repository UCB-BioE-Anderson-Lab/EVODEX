import csv
import os
from evodex.formula import calculate_exact_mass
from pipeline.version import __version__
import time
import sys
csv.field_size_limit(sys.maxsize)

BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
ERRORS_DIR = os.path.join(DATA_DIR, "errors")
EVODEX_DATA_DIR = os.path.join(BASE_DIR, "evodex", "data")

# Inputs (from Phase 3d)
EVODEX_F_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_f_phase3d_final.csv")
EVODEX_P_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_p_phase3d_final.csv")

# Outputs (pipeline)
EVODEX_M = os.path.join(PROCESSED_DIR, "EVODEX-M_unique_masses.csv")
EVODEX_M_SUBSET = os.path.join(PROCESSED_DIR, "EVODEX-M_mass_spec_subset.csv")

# Phase 5: Mass Spec Subset for EVODEX.2
# This phase generates EVODEX_M and EVODEX_M_SUBSET tables.
# EVODEX_M contains exact mass differences for all formula differences in EVODEX_F.
# EVODEX_M_SUBSET filters these to reaction patterns (EVODEX_P) compatible with single-species
# mass spec observations (i.e., at least one side of the reaction is a single fragment).
# EVODEX_M IDs are now EVODEX.2-M* using __version__ from pipeline.version.

def ensure_directories():
    for path in [PROCESSED_DIR, ERRORS_DIR, EVODEX_DATA_DIR]:
        os.makedirs(path, exist_ok=True)

def print_stats(label, value):
    print(f"{label}: {value:,}")

def parse_formula_diff_string(formula_str):
    """Parse a formula_diff string like 'C+1;H-2;N+3' into a dict.

    Returns a mapping {element: int_delta} suitable for calculate_exact_mass.
    """
    formula_str = (formula_str or "").strip()
    if not formula_str:
        return {}

    atom_diff = {}
    for term in formula_str.split(';'):
        term = term.strip()
        if not term:
            continue
        # Split element symbol (letters) from signed integer part
        i = 0
        while i < len(term) and term[i].isalpha():
            i += 1
        element = term[:i]
        delta_str = term[i:]
        if not element or not delta_str:
            raise ValueError(f"Malformed formula_diff term: {term!r}")
        delta = int(delta_str)  # e.g. '+1', '-14'
        atom_diff[element] = atom_diff.get(element, 0) + delta
    return atom_diff

def main():
    start_time = time.time()
    print("Phase 5 mass subset generation started...")

    ensure_directories()

    total_f_rows = 0
    total_m_rows = 0
    total_subset_rows = 0
    errors_f = 0
    errors_p = 0

    print("Starting Phase 5: Mass spec subset generation...")
    
    # Step 1: generate evodex_m
    evodex_m_map = {}
    with open(EVODEX_F_PHASE3D_FINAL, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            total_f_rows += 1
            try:
                atom_diff = parse_formula_diff_string(row['formula_diff'])
                mass_diff = calculate_exact_mass(atom_diff)
                evodex_m_map[row['id']] = {
                    'mass': mass_diff,
                    'sources': set(row['sources'].split(',')),
                    # Store a repr of the atom_diff dict under the original name "formula"
                    'formula': repr(atom_diff),
                }
                total_m_rows += 1
            except Exception:
                errors_f += 1
                continue

    # Step 1b: assign EVODEX_M IDs
    sorted_m_entries = sorted(evodex_m_map.items(), key=lambda x: len(x[1]['sources']), reverse=True)
    m_id_map = {}
    for i, (orig_id, data) in enumerate(sorted_m_entries, start=1):
        assigned_id = f"EVODEX.{__version__}-M{i}"
        m_id_map[orig_id] = assigned_id

    # Write EVODEX_M with new IDs (matching original schema: formula as a stringified dict)
    with open(EVODEX_M, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'mass', 'sources', 'formula'])
        writer.writeheader()
        for orig_id, data in sorted_m_entries:
            assigned_id = m_id_map[orig_id]
            writer.writerow({
                'id': assigned_id,
                'mass': data['mass'],
                'sources': ','.join(sorted(data['sources'])),
                'formula': data['formula'],
            })

    # Step 2: filter evodex_p for mass-spec compatible ones
    valid_p_ids = set()
    with open(EVODEX_P_PHASE3D_FINAL, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            try:
                left, right = row['smirks'].split('>>')
                # If you want to include potentially intramolecular reactions, use the line below instead:
                # if len(left.split('.')) == 1 or len(right.split('.')) == 1:
                if len(left.split('.')) == 1 and len(right.split('.')) == 1:
                    valid_p_ids.add(row['id'])
            except Exception:
                errors_p += 1
                continue

    # Step 3: generate evodex_m_subset
    with open(EVODEX_M_SUBSET, 'w', newline='') as outfile:
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
    dst_m = os.path.join(EVODEX_DATA_DIR, 'EVODEX-M_unique_masses.csv')
    with open(EVODEX_M, 'r') as src_file, open(dst_m, 'w', newline='') as dst_file:
        dst_file.write(src_file.read())
    print(f"Published EVODEX-M to {dst_m}")

    # EVODEX-M_SUBSET
    dst_m_subset = os.path.join(EVODEX_DATA_DIR, 'EVODEX-M_mass_spec_subset.csv')
    with open(EVODEX_M_SUBSET, 'r') as src_file, open(dst_m_subset, 'w', newline='') as dst_file:
        dst_file.write(src_file.read())
    print(f"Published EVODEX-M_SUBSET to {dst_m_subset}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 5 mass subset generation completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()