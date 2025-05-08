

import csv
import os
import hashlib
from collections import defaultdict
from evodex.formula import calculate_formula_diff
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

            evodex_id_counter = 1
            for formula_hash, data in data_map.items():
                if len(data['sources']) >= 2:
                    evodex_id = f'EVODEX.{__version__}-F{evodex_id_counter}'
                    sources = ','.join(data['sources'])
                    writer.writerow({'id': evodex_id, 'formula': str(data['formula']), 'sources': sources})
                    evodex_id_counter += 1

    return {"total": total_count, "successes": success_count, "errors": error_count}

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    print("Grouping by formula difference (EVODEX-F)...")
    result_stats = process_formula_data(
        paths['evodex_p'],
        paths['evodex_f'],
        f"{paths['errors_dir']}formula_errors.csv"
    )

    print("EVODEX-F generation complete:", result_stats)

if __name__ == "__main__":
    main()