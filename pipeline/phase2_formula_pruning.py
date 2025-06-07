import csv
import os
from collections import defaultdict
from evodex.formula import calculate_formula_diff, calculate_exact_mass
from evodex.splitting import split_reaction
from evodex.utils import reaction_hash
from pipeline.config import load_paths
from pipeline.version import __version__
import hashlib
import time

# Phase 2: Formula-Based Pruning
# This script takes the deduplicated EVODEX-R reactions and computes formula differences to derive EVODEX-F entries.
# It then filters out all EVODEX-F entries that are supported by fewer than 5 unique EVODEX-P reactions.
# Only EVODEX-F and their corresponding EVODEX-P entries that meet the support threshold are written to output.
# This pruning stage reduces noise and focuses further mining on well-supported reaction patterns.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def write_row(writer, row_data):
    writer.writerow(row_data)

def handle_error(row, e, fieldnames, error_file_path):
    error_row = {key: row[key] for key in fieldnames if key in row}
    error_row['error_message'] = str(e)
    with open(error_file_path, 'a', newline='') as errfile:
        err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
        err_writer.writerow(error_row)

def process_reaction_data(input_csv, output_csv, error_csv, process_function, additional_fields=[]):
    error_count = 0
    success_count = 0
    total_count = 0
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + additional_fields
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
            err_writer.writeheader()
        for row in reader:
            total_count += 1
            try:
                result = process_function(row)
                smirks = result["smirks"]
                if smirks.startswith(">>") or smirks.endswith(">>"):
                    continue
                writer.writerow({**row, **result})
                success_count += 1
            except Exception as e:
                error_count += 1
                handle_error(row, e, fieldnames, error_csv)
                
    return {"total":total_count, "successes":success_count, "errors":error_count}

def consolidate_reactions(input_file, output_file, prefix):
    hash_map = defaultdict(list)
    data_map = defaultdict(lambda: {'smirks': defaultdict(int), 'sources': []})
    
    with open(input_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            try:
                smirks = row['smirks']
                if smirks:  # Ensure smirks is not empty
                    rxn_hash = reaction_hash(smirks)
                    hash_map[rxn_hash].append(row['id'])
                    data_map[rxn_hash]['smirks'][smirks] += 1
                    data_map[rxn_hash]['sources'].append(row['id'])
            except Exception as e:
                pass

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        evodex_id_counter = 1
        for rxn_hash, data in data_map.items():
            if len(data['sources']) >= 2:  # Only include operators observed twice or more
                evodex_id = f'{prefix}{evodex_id_counter}'
                most_common_smirks = max(data['smirks'], key=data['smirks'].get)
                sources = ','.join(data['sources'])
                writer.writerow({'id': evodex_id, 'smirks': most_common_smirks, 'sources': sources})
                evodex_id_counter += 1

def process_formula_data(input_csv, output_csv, error_csv):
    error_count = 0
    success_count = 0
    total_count = 0
    hash_map = defaultdict(list)
    data_map = defaultdict(lambda: {'formula': None, 'sources': []})
    
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
            err_writer.writeheader()
            writer = csv.DictWriter(outfile, fieldnames=['id', 'formula', 'sources'])
            writer.writeheader()

            for row in reader:
                total_count += 1
                if total_count % 1000 == 0:
                    print(f"Processed {total_count} entries...")
                try:
                    formula_diff = calculate_formula_diff(row['smirks'])
                    formula_frozenset = frozenset(formula_diff.items())
                    formula_hash = hashlib.sha256(str(formula_frozenset).encode()).hexdigest()
                    hash_map[formula_hash].append(row['id'])
                    data_map[formula_hash]['formula'] = formula_diff
                    data_map[formula_hash]['sources'].append(row['id'])
                    success_count += 1
                except Exception as e:
                    error_count += 1
                    handle_error(row, e, fieldnames, error_csv)
            
            for formula_hash, data in data_map.items():
                writer.writerow({'id': formula_hash, 'formula': str(data['formula']), 'sources': ','.join(data['sources'])})

    return {"total": total_count, "successes": success_count, "errors": error_count}


def process_split_reactions(input_csv, output_csv, error_csv):
    error_count = 0
    success_count = 0
    total_count = 0

    hash_to_ids = defaultdict(set)
    hash_to_example_smirks = {}
    errors = []

    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=['id', 'smirks', 'sources'])
        writer.writeheader()
        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
            err_writer.writeheader()
        for row in reader:
            total_count += 1
            if total_count % 100 == 0:
                print(f"Processed {total_count} entries...")
            try:
                rxn_idx = row['id']
                smirks = row['smirks']
                split_reactions = split_reaction(smirks)
                for reaction in split_reactions:
                    reaction_hash_value = reaction_hash(reaction)
                    hash_to_ids[reaction_hash_value].add(rxn_idx)
                    if reaction_hash_value not in hash_to_example_smirks:
                        hash_to_example_smirks[reaction_hash_value] = reaction
                success_count += 1
            except Exception as e:
                error_count += 1
                handle_error(row, e, fieldnames, error_csv)
                errors.append((row['id'], str(e)))

        for reaction_hash_value, id_set in hash_to_ids.items():
            example_smirks = hash_to_example_smirks[reaction_hash_value]
            sources = ','.join(id_set)  # Ensure sources are stored as a comma-separated string
            writer.writerow({'id': reaction_hash_value, 'smirks': example_smirks, 'sources': sources})

    return {"total":total_count, "successes":success_count, "errors":error_count}


def main():
    start_time = time.time()
    print("Phase 2 formula pruning started...")
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # Step 1: Process EVODEX-R into preliminary EVODEX-P with hash as id
    process_split_reactions(paths['evodex_r_preliminary'], paths['evodex_p_preliminary'], f"{paths['errors_dir']}split_reactions_errors.csv")

    # Step 2: Process preliminary EVODEX-P into preliminary EVODEX-F with formula hash as id
    process_formula_data(paths['evodex_p_preliminary'], paths['evodex_f_preliminary'], f"{paths['errors_dir']}formula_errors.csv")

    # Step 3: Prune EVODEX-F to retain only those with >=5 sources
    retained_formula_hashes = set()
    formula_source_map = {}
    with open(paths['evodex_f_preliminary'], 'r') as f_file:
        reader = csv.DictReader(f_file)
        for row in reader:
            sources = row['sources'].split(',')
            if len(sources) >= 5:
                retained_formula_hashes.add(row['id'])
                formula_source_map[row['id']] = set(sources)

    with open(paths['evodex_f_preliminary'], 'r') as infile, open(paths['evodex_f_filtered'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if row['id'] in retained_formula_hashes:
                writer.writerow(row)

    # Step 4: Prune EVODEX-P based on retained EVODEX-F sources
    with open(paths['evodex_p_preliminary'], 'r') as infile, open(paths['evodex_p_filtered'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            for formula_sources in formula_source_map.values():
                if row['id'] in formula_sources:
                    writer.writerow(row)
                    break

    print("Phase 2 formula pruning complete.")
    print(f"Retained {len(retained_formula_hashes)} formula groups with ≥5 sources.")

    # Generate Phase 2 mining report
    print("Generating Phase 2 mining report...")

    import statistics

    # Gather EVODEX-P preliminary stats
    p_counts = []
    with open(paths['evodex_p_preliminary'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            sources = row['sources'].split(',')
            p_counts.append(len(sources))
    num_p_total = len(p_counts)

    # Gather EVODEX-F preliminary stats
    f_counts = []
    with open(paths['evodex_f_preliminary'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            sources = row['sources'].split(',')
            f_counts.append(len(sources))
    num_f_total = len(f_counts)

    # Gather retained EVODEX-F stats
    retained_f_counts = []
    with open(paths['evodex_f_filtered'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            sources = row['sources'].split(',')
            retained_f_counts.append(len(sources))
    num_f_retained = len(retained_f_counts)

    # Gather retained EVODEX-P stats
    retained_p_counts = []
    with open(paths['evodex_p_filtered'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            sources = row['sources'].split(',')
            retained_p_counts.append(len(sources))
    num_p_retained = len(retained_p_counts)

    # Compression rates
    f_compression = 100 * (num_f_total - num_f_retained) / num_f_total if num_f_total else 0
    p_compression = 100 * (num_p_total - num_p_retained) / num_p_total if num_p_total else 0

    report_lines = [
        f"EVODEX Phase 2 Formula Pruning Report (version {__version__})",
        "=============================================================",
        "",
        f"EVODEX-P Preliminary:",
        f"Total EVODEX-P generated: {num_p_total}",
        f"Min sources per EVODEX-P: {min(p_counts) if p_counts else 0}",
        f"Max sources per EVODEX-P: {max(p_counts) if p_counts else 0}",
        f"Mean sources per EVODEX-P: {statistics.mean(p_counts) if p_counts else 0:.2f}",
        f"Median sources per EVODEX-P: {statistics.median(p_counts) if p_counts else 0}",
        "",
        f"EVODEX-F Preliminary:",
        f"Total EVODEX-F generated: {num_f_total}",
        f"Min sources per EVODEX-F: {min(f_counts) if f_counts else 0}",
        f"Max sources per EVODEX-F: {max(f_counts) if f_counts else 0}",
        f"Mean sources per EVODEX-F: {statistics.mean(f_counts) if f_counts else 0:.2f}",
        f"Median sources per EVODEX-F: {statistics.median(f_counts) if f_counts else 0}",
        "",
        f"Pruning (threshold ≥5 sources):",
        f"EVODEX-F retained: {num_f_retained} of {num_f_total} ({100 - f_compression:.2f}% retained, {f_compression:.2f}% pruned)",
        f"EVODEX-P retained: {num_p_retained} of {num_p_total} ({100 - p_compression:.2f}% retained, {p_compression:.2f}% pruned)",
        "",
    ]

    # Print to console
    for line in report_lines:
        print(line)

    # Write report to file
    report_path = os.path.join(paths['errors_dir'], 'Phase2_evodex_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write('\n'.join(report_lines))

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 2 formula pruning completed in {elapsed_time:.2f} seconds.")


if __name__ == "__main__":
    main()
