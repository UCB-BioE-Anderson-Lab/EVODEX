import csv
import os
import logging
from collections import defaultdict
from evodex.decofactor import remove_cofactors
from evodex.operators import extract_operator
from evodex.formula import calculate_formula_diff, calculate_exact_mass
from evodex.utils import reaction_hash
from pipeline.config import load_paths

def setup_logging():
    logging.basicConfig(filename='mining_block.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def write_row(writer, row_data):
    writer.writerow(row_data)

def handle_error(row, e, fieldnames, error_file):
    error_row = {key: row[key] for key in fieldnames if key in row}
    error_row['error_message'] = str(e)
    logging.error(f"Error processing ID {row.get('id', 'Unknown ID')}: {e}")
    with open(error_file, 'a', newline='') as errfile:
        err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
        err_writer.writerow(error_row)

def process_reaction_data(input_csv, output_csv, error_csv, process_function, additional_fields=[]):
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + additional_fields
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
            err_writer.writeheader()
        for row in reader:
            try:
                result = process_function(row)
                writer.writerow({**row, **result})
            except Exception as e:
                handle_error(row, e, fieldnames, error_csv)

def consolidate_reactions(input_file, output_file, prefix):
    """Consolidate similar reactions based on reaction_hash."""
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
                logging.error(f"Error processing row {row['id']}: {e}")

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        evodex_id_counter = 1
        for rxn_hash, data in data_map.items():
            evodex_id = f'{prefix}{evodex_id_counter}'
            most_common_smirks = max(data['smirks'], key=data['smirks'].get)
            sources = ','.join(data['sources'])
            writer.writerow({'id': evodex_id, 'smirks': most_common_smirks, 'sources': sources})
            evodex_id_counter += 1

def process_formula_data(input_csv, output_csv, error_csv):
    """Process and consolidate formula data."""
    hash_map = defaultdict(list)
    data_map = defaultdict(lambda: {'formula': None, 'sources': []})
    
    with open(input_csv, 'r') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
            err_writer.writeheader()
        for row in reader:
            try:
                formula_diff = calculate_formula_diff(row['smirks'])
                formula_hash = hash(tuple(sorted(formula_diff.items())))
                hash_map[formula_hash].append(row['id'])
                data_map[formula_hash]['formula'] = formula_diff
                data_map[formula_hash]['sources'].append(row['id'])
            except Exception as e:
                handle_error(row, e, fieldnames, error_csv)

    with open(output_csv, 'w', newline='') as outfile:
        fieldnames = ['id', 'formula', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        evodex_id_counter = 1
        for formula_hash, data in data_map.items():
            evodex_id = f'EVODEX-F{evodex_id_counter}'
            sources = ','.join(data['sources'])
            writer.writerow({'id': evodex_id, 'formula': str(data['formula']), 'sources': sources})
            evodex_id_counter += 1

def process_mass_data(formula_csv, output_csv, error_csv):
    """Process mass data based on formula data."""
    with open(formula_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames + ['mass'])
        writer.writeheader()
        evodex_id_counter = 1
        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
            err_writer.writeheader()
        for row in reader:
            try:
                formula_diff = eval(row['formula'])
                mass_diff = calculate_exact_mass(formula_diff)
                sources = row['sources']
                evodex_id = f'EVODEX-M{evodex_id_counter}'
                writer.writerow({'id': evodex_id, 'mass': mass_diff, 'sources': sources})
                evodex_id_counter += 1
            except Exception as e:
                handle_error(row, e, fieldnames, error_csv)

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)
    setup_logging()

    def process_partial_reactions(row):
        return {'smirks': remove_cofactors(row['smirks'])}

    def process_operators(row, params):
        return {'smirks': extract_operator(row['smirks'], **params)}

    # Process partial reactions and consolidate
    process_reaction_data(paths['evodex_r'], paths['evodex_p'], f"{paths['errors_dir']}partial_reactions_errors.csv", process_partial_reactions, ['smirks'])
    print("Partial reactions processed.")
    consolidate_reactions(paths['evodex_p'], paths['evodex_p'], 'EVODEX-P')
    print("Partial reactions consolidated.")

    extract_params_map = {
        'evodex_e': {
            'include_stereochemistry': True,
            'include_sigma': True,
            'include_pi': True,
            'include_unmapped_hydrogens': True,
            'include_unmapped_heavy_atoms': True,
            'include_static_hydrogens': False
        },
        'evodex_c': {
            'include_stereochemistry': True,
            'include_sigma': False,
            'include_pi': False,
            'include_unmapped_hydrogens': True,
            'include_unmapped_heavy_atoms': True,
            'include_static_hydrogens': False
        },
        'evodex_n': {
            'include_stereochemistry': True,
            'include_sigma': True,
            'include_pi': False,
            'include_unmapped_hydrogens': True,
            'include_unmapped_heavy_atoms': True,
            'include_static_hydrogens': False
        },
        'evodex_em': {
            'include_stereochemistry': False,
            'include_sigma': True,
            'include_pi': True,
            'include_unmapped_hydrogens': False,
            'include_unmapped_heavy_atoms': False,
            'include_static_hydrogens': False
        },
        'evodex_cm': {
            'include_stereochemistry': False,
            'include_sigma': False,
            'include_pi': False,
            'include_unmapped_hydrogens': False,
            'include_unmapped_heavy_atoms': False,
            'include_static_hydrogens': False
        },
        'evodex_nm': {
            'include_stereochemistry': False,
            'include_sigma': True,
            'include_pi': False,
            'include_unmapped_hydrogens': False,
            'include_unmapped_heavy_atoms': False,
            'include_static_hydrogens': False
        }
    }

    # Process operator data and consolidate
    for key, params in extract_params_map.items():
        prefix = f'EVODEX-{key.split("_")[1].upper()}' if 'cm' not in key and 'em' not in key and 'nm' not in key else f'EVODEX-{key.split("_")[1].capitalize()}'
        error_log = f"{paths['errors_dir']}{key}_errors.csv"
        process_reaction_data(paths['evodex_p'], paths[key], error_log, lambda row: process_operators(row, params), ['smirks'])
        print(f"Operator data processed and saved to {paths[key]}")
        consolidate_reactions(paths[key], paths[key], prefix)
        print(f"Operator data consolidated and saved to {paths[key]}")

    # Process and consolidate formula data
    process_formula_data(paths['evodex_p'], paths['evodex_f'], f"{paths['errors_dir']}formula_errors.csv")
    print(f"Formula data processed and saved to {paths['evodex_f']}")

    # Process mass data based on formula data
    process_mass_data(paths['evodex_f'], paths['evodex_m'], f"{paths['errors_dir']}mass_errors.csv")
    print(f"Mass data processed and saved to {paths['evodex_m']}")

if __name__ == "__main__":
    main()
