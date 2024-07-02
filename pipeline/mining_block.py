import csv
import os
import logging
from evodex.decofactor import remove_cofactors
from evodex.operators import extract_operator
from evodex.formula import calculate_formula_diff, calculate_exact_mass
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

def handle_error(row, e):
    row['error'] = str(e)
    logging.error(f"Error processing ID {row.get('id', 'Unknown ID')}: {e}")
    return row

def process_reaction_data(input_csv, output_csv, process_function):
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['error'] if 'error' not in reader.fieldnames else reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            if row.get('error'):
                writer.writerow(row)  # Propagate rows with existing errors
                continue
            try:
                result = process_function(row)
                writer.writerow({**row, **result, 'error': ''})
            except Exception as e:
                writer.writerow(handle_error(row, e))

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)
    setup_logging()

    def process_partial_reactions(row):
        return {'smirks': remove_cofactors(row['smirks'])}

    def process_operators(row, params):
        return {'smirks': extract_operator(row['smirks'], **params)}

    def process_formula_mass(row):
        formula_diff = calculate_formula_diff(row['smirks'])
        mass_diff = calculate_exact_mass(formula_diff)
        return {'formula': str(formula_diff), 'mass': mass_diff}

    process_reaction_data(paths['evodex_r'], paths['evodex_p'], process_partial_reactions)
    print("Partial reactions processed.")

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

    for key, params in extract_params_map.items():
        process_reaction_data(paths['evodex_p'], paths[key], lambda row: process_operators(row, params))
        print(f"Operator data processed and saved to {paths[key]}")

    process_reaction_data(paths['evodex_p'], paths['evodex_f'], process_formula_mass)
    print(f"Formula data processed and saved to {paths['evodex_f']}")
    
    process_reaction_data(paths['evodex_p'], paths['evodex_m'], process_formula_mass)
    print(f"Mass data processed and saved to {paths['evodex_m']}")

if __name__ == "__main__":
    main()
