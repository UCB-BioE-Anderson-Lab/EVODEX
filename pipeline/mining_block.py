import csv
import os
import json
from collections import defaultdict
from evodex.decofactor import remove_cofactors
from evodex.operators import extract_operator
from evodex.formula import calculate_formula_diff, calculate_exact_mass
from evodex.splitting import split_reaction
from evodex.utils import reaction_hash
from pipeline.config import load_paths

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
            evodex_id = f'{prefix}{evodex_id_counter}'
            most_common_smirks = max(data['smirks'], key=data['smirks'].get)
            sources = ','.join(data['sources'])
            writer.writerow({'id': evodex_id, 'smirks': most_common_smirks, 'sources': sources})
            evodex_id_counter += 1

def process_formula_data(input_csv, output_csv, error_csv):
    hash_map = defaultdict(list)
    data_map = defaultdict(lambda: {'formula': None, 'sources': []})
    
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
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

def process_split_reactions(input_csv, output_csv, error_csv):
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
            except Exception as e:
                handle_error(row, e, fieldnames, error_csv)
                errors.append((row['id'], str(e)))

        evodex_id_counter = 1
        for reaction_hash_value, id_set in hash_to_ids.items():
            example_smirks = hash_to_example_smirks[reaction_hash_value]
            sources = json.dumps(list(id_set))
            new_id = f'EVODEX-P{evodex_id_counter}'
            writer.writerow({'id': new_id, 'smirks': example_smirks, 'sources': sources})
            evodex_id_counter += 1

import csv

def generate_synthesis_subset(input_csv, evodex_p_csv, evodex_e_csv, output_csv, error_csv):
    evodex_p_map = {}
    evodex_e_map = {}

    with open(evodex_p_csv, 'r') as p_file, open(evodex_e_csv, 'r') as e_file:
        p_reader = csv.DictReader(p_file)
        e_reader = csv.DictReader(e_file)

        # Map reaction hashes to EVODEX-P IDs
        for row in p_reader:
            reaction_hash_value = reaction_hash(row['smirks'])
            evodex_p_map[reaction_hash_value] = row['id']

        # Map EVODEX-P IDs to EVODEX-E IDs using the sources column
        for row in e_reader:
            sources = row['sources'].strip("[]").split(',')
            for source in sources:
                source = source.strip().strip('"').strip("'")
                if source in evodex_e_map:
                    evodex_e_map[source].append(row['id'])
                else:
                    evodex_e_map[source] = [row['id']]

    evodex_e_subset = set()
    evodex_e_sources = {}

    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=['id', 'sources'])
        writer.writeheader()
        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
            err_writer.writeheader()
        for row in reader:
            try:
                smirks = row['smirks']
                partial_reaction = remove_cofactors(smirks)
                reaction_hash_value = reaction_hash(partial_reaction)
                print(f"Partial reaction: {partial_reaction}, Hash: {reaction_hash_value}")
                evodex_p_id = evodex_p_map.get(reaction_hash_value)
                print(f"EVODEX-P ID: {evodex_p_id}")
                if evodex_p_id:
                    evodex_e_ids = evodex_e_map.get(evodex_p_id)
                    print(f"EVODEX-E IDs: {evodex_e_ids}")
                    if evodex_e_ids:
                        for evodex_e_id in evodex_e_ids:
                            evodex_e_subset.add(evodex_e_id)
                            if evodex_e_id in evodex_e_sources:
                                evodex_e_sources[evodex_e_id].add(evodex_p_id)
                            else:
                                evodex_e_sources[evodex_e_id] = {evodex_p_id}
            except Exception as e:
                handle_error(row, e, reader.fieldnames, error_csv)

        for evodex_e_id in sorted(evodex_e_subset):
            sources = ','.join(sorted(evodex_e_sources[evodex_e_id]))
            writer.writerow({'id': evodex_e_id, 'sources': sources})

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    def process_partial_reactions(row):
        return {'smirks': remove_cofactors(row['smirks'])}

    def process_operators(row, params):
        return {'smirks': extract_operator(row['smirks'], **params)}

    # Process split reactions and consolidate
    process_split_reactions(paths['evodex_r'], paths['evodex_p'], f"{paths['errors_dir']}split_reactions_errors.csv")

    # Generate EVODEX-E subset using decofactor
    generate_synthesis_subset(paths['evodex_r'], paths['evodex_p'], paths['evodex_e'], paths['evodex_e_synthesis'], f"{paths['errors_dir']}decofactor_errors.csv")

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
        consolidate_reactions(paths[key], paths[key], prefix)

    # Process and consolidate formula data
    process_formula_data(paths['evodex_p'], paths['evodex_f'], f"{paths['errors_dir']}formula_errors.csv")

    # Process mass data based on formula data
    process_mass_data(paths['evodex_f'], paths['evodex_m'], f"{paths['errors_dir']}mass_errors.csv")

if __name__ == "__main__":
    main()
