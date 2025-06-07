import csv
import os
import time
from collections import defaultdict
from evodex.operators import extract_operator
from evodex.utils import reaction_hash
from pipeline.config import load_paths
from pipeline.version import __version__

# Phase 3a: EVODEX-E Mining
#
# This phase extracts EVODEX-E operators from the filtered EVODEX-P set.
# Operators supported by more than one source are retained and saved for further processing.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    start_time = time.time()
    print("Phase 3a EVODEX-E mining started...")
    print("Starting Phase 3a EVODEX-E mining...")
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

        # Extract operators from EVODEX-P and build operator_map
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

    # Retain operators supported by more than 1 EVODEX-P source
    retained_ops = {k: v for k, v in operator_map.items() if len(v['sources']) > 1}

    # Write retained EVODEX-E operators to preliminary file
    with open(paths['evodex_e_phase3a_preliminary'], 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for op_hash, data in retained_ops.items():
            writer.writerow({
                'id': op_hash,
                'smirks': data['smirks'],
                'sources': ','.join(sorted(data['sources']))
            })

    # Collect statistics on the number of sources per retained operator
    import statistics
    num_sources_list = [len(data['sources']) for data in retained_ops.values()]
    if num_sources_list:
        min_sources = min(num_sources_list)
        max_sources = max(num_sources_list)
        mean_sources = statistics.mean(num_sources_list)
        median_sources = statistics.median(num_sources_list)
    else:
        min_sources = max_sources = mean_sources = median_sources = 0

    # Print statistics to the console
    print("Statistics for retained EVODEX-E operators:")
    print(f"  Number of operators: {len(retained_ops)}")
    print(f"  Sources per operator: min={min_sources}, max={max_sources}, mean={mean_sources:.2f}, median={median_sources}")

    # Save a simple report to phase3a_evodex_report.txt in the errors_dir path
    report_path = os.path.join(paths['errors_dir'], 'phase3a_evodex_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write("EVODEX-E Operator Retention Report\n")
        report_file.write(f"Number of operators: {len(retained_ops)}\n")
        report_file.write(f"Sources per operator:\n")
        report_file.write(f"  min: {min_sources}\n")
        report_file.write(f"  max: {max_sources}\n")
        report_file.write(f"  mean: {mean_sources:.2f}\n")
        report_file.write(f"  median: {median_sources}\n")

    print(f"Phase 3a EVODEX-E mining complete. Retained {len(retained_ops)} operators with >1 source.")

    # === Prune EVODEX-P to only retained P entries ===
    # Build set of all surviving P IDs from retained_ops['sources']
    surviving_p_ids = set()
    for data in retained_ops.values():
        surviving_p_ids.update(data['sources'])

    # Load EVODEX-P filtered entries again
    p_rows = []
    with open(paths['evodex_p_filtered'], 'r') as pfile:
        reader = csv.DictReader(pfile)
        for row in reader:
            p_rows.append(row)

    # Write new CSV to paths['evodex_p_phase3a_retained'] with only P rows whose ID is in the surviving set
    retained_p_rows = [row for row in p_rows if row['id'] in surviving_p_ids]
    with open(paths['evodex_p_phase3a_retained'], 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=p_rows[0].keys())
        writer.writeheader()
        writer.writerows(retained_p_rows)

    # Print pruning summary
    print(f"Pruned EVODEX-P: {len(retained_p_rows)} retained of {len(p_rows)} original entries.")

    # After processing all operators, print number of extraction errors recorded
    try:
        with open(error_log_path, 'r') as errfile:
            error_count = sum(1 for _ in errfile) - 1  # subtract header
    except Exception:
        error_count = 0
    print(f"Number of extraction errors recorded: {error_count}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 3a EVODEX-E mining completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()