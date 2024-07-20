import csv
import os
from collections import defaultdict
from evodex.filter import validate_smiles
from evodex.astatine import hydrogen_to_astatine_reaction
from evodex.mapping import map_atoms
from evodex.utils import reaction_hash
from pipeline.config import load_paths

def ensure_directories(paths: dict):
    """Ensure that all necessary directories exist."""
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def write_row(writer, data):
    """Write a row to the CSV file."""
    writer.writerow(data)

def process_raw_data(input_file, output_file, ec_representation):
    """Process the initial raw data file."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'smirks', 'sources', 'error', 'ec_num']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            ec_num = row.get('ec_num', '')
            if ec_num:
                ec_representation[ec_num] = False
            try:
                if validate_smiles(row['mapped']):
                    new_row = {'id': row['rxn_idx'], 'smirks': row['mapped'], 'sources': row.get('rxn_idx', ''), 'error': '', 'ec_num': ec_num}
                else:
                    raise ValueError("SMILES validation failed")
            except Exception as e:
                new_row = {'id': row['rxn_idx'], 'smirks': '', 'sources': '', 'error': str(e), 'ec_num': ec_num}
            write_row(writer, new_row)


def process_data(input_file, output_file, transformation_function):
    """General function to process data with transformation."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'smirks', 'sources', 'error', 'ec_num']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            if row['error']:  # Skip processing if there's an existing error
                write_row(writer, row)
                continue
            try:
                new_data = transformation_function(row)
                write_row(writer, new_data)
            except Exception as e:
                row['error'] = str(e)
                write_row(writer, row)

def consolidate_reactions(input_file, output_file, ec_representation):
    """Consolidate similar reactions based on reaction_hash and retain up to 2 unique instances of each ec_num."""
    ec_num_map = defaultdict(lambda: defaultdict(dict))

    with open(input_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            try:
                smirks = row['smirks']
                ec_num = row.get('ec_num', '')
                if smirks and ec_num and len(ec_num_map[ec_num]) < 2:  # Ensure smirks and ec_num are not empty and limit to 2 instances
                    rxn_hash = reaction_hash(smirks)
                    if rxn_hash not in ec_num_map[ec_num]:  # Ensure unique reaction hash
                        ec_num_map[ec_num][rxn_hash] = row
            except Exception as e:
                print(f"Error processing row: {row} - {e}")

    sorted_ec_nums = sorted(ec_num_map.keys(), key=lambda ec: list(map(int, ec.split('.'))))

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources', 'ec_num', 'error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        evodex_id_counter = 1
        for ec_num in sorted_ec_nums:
            for rxn_hash, data in ec_num_map[ec_num].items():
                evodex_id = f'EVODEX-R{evodex_id_counter}'
                data['id'] = evodex_id
                sources = data['sources']
                writer.writerow({
                    'id': evodex_id,
                    'smirks': data['smirks'],
                    'sources': sources,
                    'ec_num': data['ec_num'],
                    'error': data['error']
                })
                ec_representation[ec_num] = True
                evodex_id_counter += 1
    print(f"Total reactions written: {evodex_id_counter - 1}")

def write_ec_representation(ec_representation, output_file):
    """Write the EC number representation dictionary to a CSV file."""
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['ec_num', 'represented'])
        for ec_num, represented in ec_representation.items():
            writer.writerow([ec_num, represented])

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    ec_representation = {}

    # Process initial raw data
    process_raw_data(paths['raw_data'], paths['filtered_data'], ec_representation)

    # Subsequent processing steps
    process_data(paths['filtered_data'], paths['astatine_data'], lambda row: {
        'id': row['id'],
        'smirks': hydrogen_to_astatine_reaction(row['smirks']),
        'sources': row['sources'],
        'ec_num': row['ec_num'],
        'error': ''
    })

    process_data(paths['astatine_data'], paths['mapped_data'], lambda row: {
        'id': row['id'],
        'smirks': map_atoms(row['smirks']),
        'sources': row['sources'],
        'ec_num': row['ec_num'],
        'error': ''
    })

    # Consolidate reactions to generate EVODEX-R data with EC number constraints
    consolidate_reactions(paths['mapped_data'], paths['evodex_r'], ec_representation)

    # Write EC number representation to a CSV file
    write_ec_representation(ec_representation, os.path.join(paths['errors_dir'], 'ec_representation.csv'))

if __name__ == "__main__":
    main()
