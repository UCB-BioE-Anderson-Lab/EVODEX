import csv
import os
import logging
from evodex.filter import validate_smiles
from evodex.astatine import hydrogen_to_astatine_reaction
from evodex.mapping import map_atoms
from pipeline.config import load_paths

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(filename='data_processing.log', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

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

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)
    setup_logging()

    # Process initial raw data
    process_raw_data(paths['raw_data'], paths['filtered_data'])
    print(f"Filtered reaction data generated and saved to {paths['filtered_data']}")

    # Subsequent processing steps
    process_data(paths['filtered_data'], paths['astatine_data'], lambda row: {
        'id': row['id'],
        'smirks': hydrogen_to_astatine_reaction(row['smirks']),
        'sources': row['sources'],
        'error': ''
    })
    print(f"Astatine test data generated and saved to {paths['astatine_data']}")

    process_data(paths['astatine_data'], paths['mapped_data'], lambda row: {
        'id': row['id'],
        'smirks': map_atoms(row['smirks']),
        'sources': row['sources'],
        'error': ''
    })
    print(f"Atom-mapped test data generated and saved to {paths['mapped_data']}")

if __name__ == "__main__":
    main()
