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

def process_raw_data(input_file, output_file):
    """Process the initial raw data file."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'smirks', 'sources', 'error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            try:
                if validate_smiles(row['mapped']):
                    new_row = {'id': row['rxn_idx'], 'smirks': row['mapped'], 'sources': row.get('sources', ''), 'error': ''}
                else:
                    raise ValueError("SMILES validation failed")
            except Exception as e:
                new_row = {'id': row['rxn_idx'], 'smirks': '', 'sources': '', 'error': str(e)}
            write_row(writer, new_row)

def process_data(input_file, output_file, transformation_function):
    """General function to process data with transformation."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'smirks', 'sources', 'error']
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

def consolidate_reactions(input_file, output_file):
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
                pass

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        evodex_id_counter = 1
        for rxn_hash, data in data_map.items():
            evodex_id = f'EVODEX-R{evodex_id_counter}'
            most_common_smirks = max(data['smirks'], key=data['smirks'].get)
            sources = ','.join(data['sources'])
            writer.writerow({'id': evodex_id, 'smirks': most_common_smirks, 'sources': sources})
            evodex_id_counter += 1

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # Process initial raw data
    process_raw_data(paths['raw_data'], paths['filtered_data'])

    # Subsequent processing steps
    process_data(paths['filtered_data'], paths['astatine_data'], lambda row: {
        'id': row['id'],
        'smirks': hydrogen_to_astatine_reaction(row['smirks']),
        'sources': row['sources'],
        'error': ''
    })

    process_data(paths['astatine_data'], paths['mapped_data'], lambda row: {
        'id': row['id'],
        'smirks': map_atoms(row['smirks']),
        'sources': row['sources'],
        'error': ''
    })

    # Consolidate reactions to generate EVODEX-R data
    consolidate_reactions(paths['mapped_data'], paths['evodex_r'])

if __name__ == "__main__":
    main()
