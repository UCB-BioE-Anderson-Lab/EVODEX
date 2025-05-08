

import csv
import os
from collections import defaultdict
from evodex.mapping import mine_operator_from_mapping
from evodex.utils import hash_operator
from pipeline.config import load_paths
from pipeline.version import __version__

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def handle_error(row, e, fieldnames, error_file_path):
    error_row = {key: row.get(key, "") for key in fieldnames}
    error_row['error_message'] = str(e)
    with open(error_file_path, 'a', newline='') as errfile:
        err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
        err_writer.writerow(error_row)

def mine_eros(input_csv, output_csv, error_csv):
    error_count = 0
    success_count = 0
    total_count = 0

    seen_hashes = set()
    hash_to_sources = defaultdict(list)

    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'ero_hash', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=reader.fieldnames + ['error_message'])
            err_writer.writeheader()

            for row in reader:
                total_count += 1
                try:
                    p_id = row['id']
                    smirks = row['smirks']
                    ero_obj = mine_operator_from_mapping(smirks)
                    ero_smirks = ero_obj.as_smirks()
                    ero_hash = hash_operator(ero_obj)
                    hash_to_sources[ero_hash].append(p_id)

                    if ero_hash not in seen_hashes:
                        writer.writerow({
                            'id': f'EVODEX.{__version__}-E{len(seen_hashes)+1}',
                            'ero_hash': ero_hash,
                            'smirks': ero_smirks,
                            'sources': p_id
                        })
                        seen_hashes.add(ero_hash)

                    success_count += 1
                except Exception as e:
                    error_count += 1
                    handle_error(row, e, reader.fieldnames, error_csv)

    return {"total": total_count, "successes": success_count, "errors": error_count}

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    print("Mining initial EROs from EVODEX-P...")
    result_stats = mine_eros(
        paths['evodex_p'],
        paths['evodex_e_raw'],
        f"{paths['errors_dir']}ero_mining_errors.csv"
    )

    print("ERO mining complete:", result_stats)

if __name__ == "__main__":
    main()