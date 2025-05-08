

import csv
import os
from collections import defaultdict
from evodex.utils import reaction_hash
from evodex.splitting import split_reaction
from pipeline.config import load_paths
from pipeline.version import __version__

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def handle_error(row, e, fieldnames, error_file_path):
    error_row = {key: row[key] for key in fieldnames if key in row}
    error_row['error_message'] = str(e)
    with open(error_file_path, 'a', newline='') as errfile:
        err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
        err_writer.writerow(error_row)

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
            try:
                rxn_idx = row['id']
                smirks = row['smirks']
                split_reactions = split_reaction(smirks)
                for reaction in split_reactions:
                    reaction_hash_value = reaction_hash(reaction)
                    if reaction_hash_value not in hash_to_ids:
                        hash_to_ids[reaction_hash_value].add(rxn_idx)
                        hash_to_example_smirks[reaction_hash_value] = reaction
                    else:
                        hash_to_ids[reaction_hash_value].add(rxn_idx)
                success_count += 1
            except Exception as e:
                error_count += 1
                handle_error(row, e, fieldnames, error_csv)
                errors.append((row['id'], str(e)))

        evodex_id_counter = 1
        for reaction_hash_value, id_set in hash_to_ids.items():
            example_smirks = hash_to_example_smirks[reaction_hash_value]
            sources = ','.join(id_set)
            new_id = f'EVODEX.{__version__}-P{evodex_id_counter}'
            writer.writerow({'id': new_id, 'smirks': example_smirks, 'sources': sources})
            evodex_id_counter += 1

    return {"total": total_count, "successes": success_count, "errors": error_count}

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    print("Generating EVODEX-P...")
    result_stats = process_split_reactions(
        paths['evodex_r'],
        paths['evodex_p'],
        f"{paths['errors_dir']}split_reactions_errors.csv"
    )

    print("EVODEX-P generation complete:", result_stats)

if __name__ == "__main__":
    main()