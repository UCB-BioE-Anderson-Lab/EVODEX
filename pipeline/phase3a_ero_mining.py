import csv
import os
import time
from collections import defaultdict
from evodex.operators import extract_operator
from evodex.utils import reaction_hash
from pipeline.config import load_paths
from pipeline.version import __version__
import sys
csv.field_size_limit(sys.maxsize)

# Phase 3a: EVODEX-E Mining
# This phase extracts EVODEX-E operators from the filtered EVODEX-P set.
# Operators with >1 sources are retained. Operators with only 1 source are excluded. Up to 5 P's are kept per retained E.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    start_time = time.time()
    print("Phase 3a EVODEX-E mining started...")
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # Load EVODEX-P filtered entries and extract operators within error log context
    p_entries = []
    error_log_path = os.path.join(paths['errors_dir'], 'phase3a_ero_errors.csv')
    with open(error_log_path, 'w', newline='') as errorfile:
        err_writer = csv.DictWriter(errorfile, fieldnames=['id', 'smirks', 'error_message'])
        err_writer.writeheader()

        with open(paths['evodex_p_filtered'], 'r') as pfile:
            reader = csv.DictReader(pfile)
            for row in reader:
                p_entries.append(row)

        # Extract EVODEX-E operators and record supporting sources
        operator_map = defaultdict(lambda: {'smirks': None, 'sources': set()})
        for i, row in enumerate(p_entries, start=1):
            if i % 100 == 0:
                print(f"Processing EVODEX-P row {i} for operator extraction...")
            try:
                op_smirks = extract_operator(
                    row['smirks'],
                    include_stereochemistry=True,
                    include_sigma=True,
                    include_pi=True,
                    include_unmapped_hydrogens=True,
                    include_unmapped_heavy_atoms=True,
                    include_static_hydrogens=True
                )
                if op_smirks.startswith(">>") or op_smirks.endswith(">>"):
                    continue
                op_hash = reaction_hash(op_smirks)
                operator_map[op_hash]['smirks'] = op_smirks
                operator_map[op_hash]['sources'].add(row['id'])
            except Exception as e:
                err_writer.writerow({'id': row['id'], 'smirks': row['smirks'], 'error_message': str(e)})

    # Retain operators with >1 sources (operators with only 1 source are excluded)
    retained_ops = {k: v for k, v in operator_map.items() if len(v['sources']) > 1}

    # Statistics on retained operators
    import statistics
    num_sources_list = [len(data['sources']) for data in retained_ops.values()]
    min_sources = min(num_sources_list) if num_sources_list else 0
    max_sources = max(num_sources_list) if num_sources_list else 0
    mean_sources = statistics.mean(num_sources_list) if num_sources_list else 0
    median_sources = statistics.median(num_sources_list) if num_sources_list else 0

    # Retention breakdown (dropping operators with only 1 source, keeping up to 5 sources per retained E)
    num_ops_excluded = sum(1 for v in operator_map.values() if len(v['sources']) <= 1)
    num_ops_retained = len(retained_ops)
    num_ops_trimmed = sum(1 for v in retained_ops.values() if len(v['sources']) > 5)

    print("Operator retention breakdown:")
    print(f"  Operators with only 1 source (excluded): {num_ops_excluded}")
    print(f"  Operators retained (>1 sources): {num_ops_retained}")
    print(f"  Operators with >5 sources and thus trimmed to 5: {num_ops_trimmed}")

    # Sort retained_ops by number of sources (descending)
    sorted_retained_ops = sorted(
        retained_ops.items(),
        key=lambda item: len(item[1]['sources']),
        reverse=True
    )

    # Write retained EVODEX-E operators to preliminary file (trim sources to up to 5 per operator)
    with open(paths['evodex_e_phase3a_preliminary'], 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for op_hash, data in sorted_retained_ops:
            selected_sources = sorted(data['sources'])[:5] if data['sources'] else []
            writer.writerow({
                'id': op_hash,
                'smirks': data['smirks'],
                'sources': ','.join(selected_sources)
            })

    # Print statistics
    print("Statistics for retained EVODEX-E operators:")
    print(f"  Number of operators: {len(retained_ops)}")
    print(f"  Sources per operator: min={min_sources}, max={max_sources}, mean={mean_sources:.2f}, median={median_sources}")

    # Save report
    report_path = os.path.join(paths['errors_dir'], 'phase3a_evodex_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write("EVODEX-E Operator Retention Report\n")
        report_file.write(f"Number of operators: {len(retained_ops)}\n")
        report_file.write(f"Sources per operator:\n")
        report_file.write(f"  min: {min_sources}\n")
        report_file.write(f"  max: {max_sources}\n")
        report_file.write(f"  mean: {mean_sources:.2f}\n")
        report_file.write(f"  median: {median_sources}\n")
        report_file.write(f"\nOperator retention breakdown:\n")
        report_file.write(f"  Operators with only 1 source (excluded): {num_ops_excluded}\n")
        report_file.write(f"  Operators retained (>1 sources): {num_ops_retained}\n")
        report_file.write(f"  Operators with >5 sources and thus trimmed to 5: {num_ops_trimmed}\n")

    print(f"Phase 3a EVODEX-E mining complete. Retained {len(retained_ops)} operators with >1 sources. Up to 5 sources retained per operator.")

    # === Prune EVODEX-P to only retained P entries ===
    # Build set of surviving P IDs — up to 5 per E
    surviving_p_ids = set()
    for data in retained_ops.values():
        selected_sources = sorted(data['sources'])[:5] if data['sources'] else []
        surviving_p_ids.update(selected_sources)

    # Load EVODEX-P again
    p_rows = []
    with open(paths['evodex_p_filtered'], 'r') as pfile:
        reader = csv.DictReader(pfile)
        for row in reader:
            p_rows.append(row)

    # Filter P entries
    retained_p_rows = [row for row in p_rows if row['id'] in surviving_p_ids]

    eliminated_p_count = len(p_rows) - len(retained_p_rows)
    print(f"  EVODEX-P entries eliminated by pruning: {eliminated_p_count}")

    with open(report_path, 'a') as report_file:
        report_file.write(f"\nEVODEX-P entries eliminated by pruning: {eliminated_p_count}\n")

    # Write retained EVODEX-P
    with open(paths['evodex_p_phase3a_retained'], 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=p_rows[0].keys())
        writer.writeheader()
        writer.writerows(retained_p_rows)

    print(f"Pruned EVODEX-P: {len(retained_p_rows)} retained of {len(p_rows)} original entries.")

    # Report extraction errors
    try:
        with open(error_log_path, 'r') as errfile:
            error_count = sum(1 for _ in errfile) - 1  # subtract header
    except Exception:
        error_count = 0
    print(f"Number of extraction errors recorded: {error_count}")

    # Done
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 3a EVODEX-E mining completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()