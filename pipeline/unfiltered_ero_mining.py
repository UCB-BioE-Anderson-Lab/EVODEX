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

# Phase 3: EVODEX-E Mining
# This phase extracts all EVODEX-E operators from the filtered EVODEX-P set.
# All extracted operators and their sources are retained.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    start_time = time.time()
    print("Phase 3 EVODEX-E mining started...")
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # Load EVODEX-P filtered entries and extract operators within error log context
    p_entries = []
    error_log_path = os.path.join(paths['errors_dir'], 'phase3a_ero_errors.csv')
    with open(error_log_path, 'w', newline='') as errorfile:
        err_writer = csv.DictWriter(errorfile, fieldnames=['id', 'smirks', 'error_message'])
        err_writer.writeheader()

        with open(paths['evodex_r_preliminary'], 'r') as pfile:
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

    # Write all extracted EVODEX-E operators to output file (all sources retained)
    with open('data/processed/unfiltered_eros.csv', 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for op_hash, data in operator_map.items():
            selected_sources = sorted(data['sources']) if data['sources'] else []
            writer.writerow({
                'id': op_hash,
                'smirks': data['smirks'],
                'sources': ','.join(selected_sources)
            })

    # Statistics on extracted operators
    import statistics
    num_sources_list = [len(data['sources']) for data in operator_map.values()]
    min_sources = min(num_sources_list) if num_sources_list else 0
    max_sources = max(num_sources_list) if num_sources_list else 0
    mean_sources = statistics.mean(num_sources_list) if num_sources_list else 0
    median_sources = statistics.median(num_sources_list) if num_sources_list else 0

    print("Statistics for extracted EVODEX-E operators:")
    print(f"  Number of operators: {len(operator_map)}")
    print(f"  Sources per operator: min={min_sources}, max={max_sources}, mean={mean_sources:.2f}, median={median_sources}")

    # Save report
    report_path = os.path.join(paths['errors_dir'], 'phase3a_evodex_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write("EVODEX-E Operator Extraction Report - Phase 3 Mining\n")
        report_file.write(f"Number of operators: {len(operator_map)}\n")
        report_file.write(f"Sources per operator:\n")
        report_file.write(f"  min: {min_sources}\n")
        report_file.write(f"  max: {max_sources}\n")
        report_file.write(f"  mean: {mean_sources:.2f}\n")
        report_file.write(f"  median: {median_sources}\n")

    print(f"Phase 3 EVODEX-E mining complete. Extracted {len(operator_map)} operators with all sources retained.")

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
    print(f"Phase 3 EVODEX-E mining completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()