import csv
import os
import time
from collections import defaultdict, Counter
import statistics
from evodex.astatine import hydrogen_to_astatine_reaction
from evodex.mapping import map_atoms
from evodex.utils import reaction_hash
from pipeline.config import load_paths
from pipeline.version import __version__

# Phase 1: Data Preparation
# This script starts from the BRENDA dataset (downloaded via EnzymeMap) and reduces it to a deduplicated bag of unique reactions.
# The reactions are stripped of enzyme-specific information, retaining only substrate and product SMILES.
# Reactions come with atom maps already applied. Hydrogens are replaced with astatines, and new atom maps are assigned to the 
# astatines using the next available map index.


def ensure_directories(paths: dict):
    """Ensure that all necessary directories exist."""
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def write_row(writer, data):
    """Write a row to the CSV file."""
    writer.writerow(data)

def process_raw_data(input_file, output_file):
    """Process the initial raw data file."""
    total_raw_reactions = 0
    valid_reactions = 0

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'smirks', 'sources', 'error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            if not any(row.values()):
                continue
            total_raw_reactions += 1
            if total_raw_reactions % 100 == 0:
                print(f"[{time.strftime('%H:%M:%S')}] Processed {total_raw_reactions} raw reactions...")
            try:
                # This is a validation step to exclude any malformed data
                reaction_hash(row['mapped'])
                valid_reactions += 1

                # Write the validated data row
                new_row = {'id': row['rxn_idx'], 'smirks': row['mapped'], 'sources': row.get('rxn_idx', ''), 'error': ''}
            except Exception as e:
                # Write an erroneous data row
                new_row = {'id': row['rxn_idx'], 'smirks': '', 'sources': '', 'error': str(e)}
            write_row(writer, new_row)

    return total_raw_reactions, valid_reactions

def process_data(input_file, output_file, transformation_function, stage_name, error_log):
    """General function to process data with transformation."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'smirks', 'sources', 'error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for idx, row in enumerate(reader, start=1):
            if idx % 100 == 0:
                print(f"[{time.strftime('%H:%M:%S')}] {stage_name} processing row {idx}...")
            if row['error']:  # Skip processing if there's an existing error
                write_row(writer, row)
                continue
            try:
                new_data = transformation_function(row)
                write_row(writer, new_data)
            except Exception as e:
                import traceback
                error_msg = f"{str(e)}\n{traceback.format_exc()}"
                row['error'] = error_msg
                error_log.append({**row, 'stage': stage_name})
                write_row(writer, row)


def main():
    start_time = time.time()
    print("Phase 1 data preparation started...")

    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # Process initial raw data
    total_raw_reactions, valid_reactions = process_raw_data(paths['raw_data'], paths['filtered_data'])

    error_log = []

    # Subsequent processing steps
    process_data(paths['filtered_data'], paths['astatine_data'], lambda row: {
        'id': row['id'],
        'smirks': hydrogen_to_astatine_reaction(row['smirks']),
        'sources': row['sources'],
        'error': ''
    }, 'astatine', error_log)

    process_data(paths['astatine_data'], paths['mapped_data'], lambda row: {
        'id': row['id'],
        'smirks': map_atoms(row['smirks']),
        'sources': row['sources'],
        'error': ''
    }, 'mapping', error_log)

    # Calculate statistics
    percentage_loss = ((total_raw_reactions - valid_reactions) / total_raw_reactions) * 100

    # Generate mining report
    smirks_counts = Counter()
    with open(paths['mapped_data'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            if not row['error']:
                smirks = row['smirks']
                rxn_hash = reaction_hash(smirks)
                smirks_counts[rxn_hash] += 1

    unique_ev_r = len(smirks_counts)
    counts = list(smirks_counts.values())

    report_lines = [
        f"EVODEX-R Mining Report (version {__version__})",
        "============================================",
        f"Total raw reactions: {total_raw_reactions}",
        f"Valid reactions after sanitization: {valid_reactions}",
        f"Compression rate due to reaction deduplication: {percentage_loss:.2f}%",
        f"(i.e., only {100 - percentage_loss:.2f}% of raw reactions are unique after hashing)",
        "",
        f"Total unique EVODEX-R entries: {unique_ev_r}",
        f"Min reactions per EVODEX-R: {min(counts) if counts else 0}",
        f"Max reactions per EVODEX-R: {max(counts) if counts else 0}",
        f"Mean reactions per EVODEX-R: {statistics.mean(counts) if counts else 0:.2f}",
        f"Median reactions per EVODEX-R: {statistics.median(counts) if counts else 0}",
        "",
    ]

    for line in report_lines:
        print(line)

    report_path = os.path.join(paths['errors_dir'], 'Phase1_evodex_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write('\n'.join(report_lines))

    # Write evodex_r_raw.csv with deduplicated reactions
    unique_reactions = {}
    with open(paths['mapped_data'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            if not row['error']:
                rxn_hash = reaction_hash(row['smirks'])
                if rxn_hash not in unique_reactions:
                    unique_reactions[rxn_hash] = {
                        'id': rxn_hash,
                        'smirks': row['smirks'],
                        'sources': row['sources']
                    }

    with open(paths['evodex_r_preliminary'], 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for data in unique_reactions.values():
            writer.writerow(data)

    if error_log:
        error_file_path = os.path.join(paths['errors_dir'], 'data_preparation_errors.csv')
        fieldnames = list(error_log[0].keys())
        with open(error_file_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(error_log)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 1 data preparation completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
