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
from pipeline.version import __version__
import hashlib

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
    
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile, open(error_csv, 'w', newline='') as errfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
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
                handle_error(row, e, fieldnames, err_writer)
                
        evodex_id_counter = 1
        for formula_hash, data in data_map.items():
            if len(data['sources']) >= 2:  # Only include if observed at least twice
                evodex_id = f'EVODEX.{__version__}-F{evodex_id_counter}'
                sources = ','.join(data['sources'])
                writer.writerow({'id': evodex_id, 'formula': str(data['formula']), 'sources': sources})
                evodex_id_counter += 1

    return {"total": total_count, "successes": success_count, "errors": error_count}

def process_mass_data(formula_csv, output_csv, error_csv):
    error_count = 0
    success_count = 0
    total_count = 0

    with open(formula_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['mass']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        evodex_id_counter = 1
        with open(error_csv, 'w', newline='') as errfile:
            err_writer = csv.DictWriter(errfile, fieldnames=fieldnames + ['error_message'])
            err_writer.writeheader()
        for row in reader:
            total_count += 1
            try:
                formula_diff = eval(row['formula'])
                mass_diff = calculate_exact_mass(formula_diff)
                sources = row['sources']
                evodex_id = f'EVODEX.{__version__}-M{evodex_id_counter}'
                writer.writerow({'id': evodex_id, 'mass': mass_diff, 'sources': sources, 'formula': row['formula']})
                evodex_id_counter += 1
                success_count += 1
            except Exception as e:
                error_count += 1
                handle_error(row, e, fieldnames, error_csv)

    return {"total":total_count, "successes":success_count, "errors":error_count}

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
            sources = ','.join(id_set)  # Ensure sources are stored as a comma-separated string
            new_id = f'EVODEX.{__version__}-P{evodex_id_counter}'
            writer.writerow({'id': new_id, 'smirks': example_smirks, 'sources': sources})
            evodex_id_counter += 1

    return {"total":total_count, "successes":success_count, "errors":error_count}

def generate_synthesis_subset(input_csv, evodex_p_csv, evodex_e_csv, output_csv, error_csv):
    error_count = 0
    success_count = 0
    total_count = 0

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
                total_count += 1
                try:
                    smirks = row['smirks']
                    partial_reaction = remove_cofactors(smirks)
                    # Check if partial_reaction starts or ends with >>
                    if partial_reaction.startswith(">>") or partial_reaction.endswith(">>"):
                        continue
                    
                    reaction_hash_value = reaction_hash(partial_reaction)
                    evodex_p_id = evodex_p_map.get(reaction_hash_value)

                    if evodex_p_id:
                        evodex_e_ids = evodex_e_map.get(evodex_p_id)
                        if evodex_e_ids:
                            for evodex_e_id in evodex_e_ids:
                                evodex_e_subset.add(evodex_e_id)
                                if evodex_e_id in evodex_e_sources:
                                    evodex_e_sources[evodex_e_id].add(evodex_p_id)
                                else:
                                    evodex_e_sources[evodex_e_id] = {evodex_p_id}
                except Exception as e:
                    error_count += 1
                    handle_error(row, e, reader.fieldnames, error_csv)

        for evodex_e_id in sorted(evodex_e_subset):
            sources = evodex_e_sources[evodex_e_id]
            sources_str = ','.join(sorted(sources))
            writer.writerow({'id': evodex_e_id, 'sources': sources_str})
            success_count += 1

    return {"total":total_count, "successes":success_count, "errors":error_count}

def generate_mass_spec_subset(evodex_p_csv, evodex_m_csv, output_csv, error_csv):
    error_count = 0
    success_count = 0
    total_count = 0    
    
    valid_evodex_p_ids = set()
    evodex_p_map = {}

    with open(evodex_p_csv, 'r') as p_file:
        p_reader = csv.DictReader(p_file)
        
        for row in p_reader:
            total_count += 1
            try:
                splitted = row['smirks'].split('>>')
                substrates = splitted[0].split('.')
                products = splitted[1].split('.')

                # Check if either substrates or products have a single molecule
                if len(substrates) == 1 or len(products) == 1:
                    valid_evodex_p_ids.add(row['id'])
                    success_count += 1
                
            except Exception as e:
                error_count += 1

    # Create a map from each sources id in evodex-m df (evodex-p values) to the id field (evodex-m)
    evodex_m_to_p_map = defaultdict(set)
    evodex_m_to_mass_map = {}

    with open(evodex_m_csv, 'r') as m_file:
        m_reader = csv.DictReader(m_file)
        
        for row in m_reader:
            sources = row['sources'].split(',')
            for source in sources:
                if source in valid_evodex_p_ids:
                    evodex_m_to_p_map[row['id']].add(source)
                    evodex_m_to_mass_map[row['id']] = row['mass']

    # Write the output CSV
    with open(output_csv, 'w', newline='') as outfile:
        fieldnames = ['id', 'mass', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for evodex_m, sources in evodex_m_to_p_map.items():
            sources_str = ','.join(sources)
            mass = evodex_m_to_mass_map[evodex_m]
            writer.writerow({'id': evodex_m, 'mass': mass, 'sources': sources_str})

    return {"total":total_count, "successes":success_count, "errors":error_count}

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    def process_partial_reactions(row):
        return {'smirks': remove_cofactors(row['smirks'])}

    def process_operators(row, params):
        return {'smirks': extract_operator(row['smirks'], **params)}

    result_stats = {}

    # Process split reactions and consolidate
    result_stats["evodex_p"] = process_split_reactions(paths['evodex_r'], paths['evodex_p'], f"{paths['errors_dir']}split_reactions_errors.csv")

    extract_params_map = {
        'evodex_e': {
            'include_stereochemistry': True,
            'include_sigma': True,
            'include_pi': True,
            'include_unmapped_hydrogens': True,
            'include_unmapped_heavy_atoms': True,
            'include_static_hydrogens': True
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
            'include_static_hydrogens': True
        },
        'evodex_em': {
            'include_stereochemistry': False,
            'include_sigma': True,
            'include_pi': True,
            'include_unmapped_hydrogens': False,
            'include_unmapped_heavy_atoms': False,
            'include_static_hydrogens': True
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
            'include_static_hydrogens': True
        }
    }

    # Process operator data and consolidate
    for key, params in extract_params_map.items():
        prefix = f'EVODEX.{__version__}-{key.split("_")[1].upper()}' if 'cm' not in key and 'em' not in key and 'nm' not in key else f'EVODEX.{__version__}-{key.split("_")[1].capitalize()}'
        error_log = f"{paths['errors_dir']}{key}_errors.csv"
        result = process_reaction_data(paths['evodex_p'], paths[key], error_log, lambda row: process_operators(row, params), ['smirks'])
        result_stats[key] = result
        consolidate_reactions(paths[key], paths[key], prefix)

    # Process and consolidate formula data
    result_stats["evodex_f"] = process_formula_data(paths['evodex_p'], paths['evodex_f'], f"{paths['errors_dir']}formula_errors.csv")

    # Process mass data based on formula data
    result_stats["evodex_m"] = process_mass_data(paths['evodex_f'], paths['evodex_m'], f"{paths['errors_dir']}mass_errors.csv")

    # Generate EVODEX-E subset using decofactor
    result_stats["synthesis"] = generate_synthesis_subset(paths['evodex_r'], paths['evodex_p'], paths['evodex_e'], paths['evodex_e_synthesis'], f"{paths['errors_dir']}decofactor_errors.csv")

    # Generate mass spec subset
    result_stats["mass_spec"] = generate_mass_spec_subset(paths['evodex_p'], paths['evodex_m'], paths['evodex_m_subset'], f"{paths['errors_dir']}mass_spec_subset_errors.csv")

    def write_result_stats_to_csv(stats, output_file):
        fieldnames = ['process', 'total', 'successes', 'errors']
        
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for process, data in stats.items():
                row = {'process': process, **data}
                writer.writerow(row)

    # Usage
    output_file = 'data/errors/result_stats.csv'
    write_result_stats_to_csv(result_stats, output_file)

    print(result_stats)


if __name__ == "__main__":
    main()
