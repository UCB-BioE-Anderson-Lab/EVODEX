import csv
import os
from evodex.filter import validate_smiles
from evodex.astatine import hydrogen_to_astatine_reaction
from evodex.mapping import map_atoms
from pipeline.config import load_paths

def generate_filtered_reaction_data(input_file: str, output_file: str):
    """Generate a CSV file with validated reactions from the input file."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for row in reader:
            smiles = row['mapped']  # Use the correct column name for reaction SMILES
            if validate_smiles(smiles):
                writer.writerow(row)

def generate_astatine_data(input_file: str, output_file: str):
    """Generate a CSV file with astatine-containing reactions from hydrogen-containing reactions."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['rxn_idx', 'astatine_mapped']  # Assuming 'rxn_idx' is the ID column
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for row in reader:
            reaction_id = row['rxn_idx']  # Use the correct column name for reaction ID
            hydrogen_smiles = row['mapped']  # Use the correct column name for hydrogen SMILES
            astatine_smiles = hydrogen_to_astatine_reaction(hydrogen_smiles)
            new_row = {'rxn_idx': reaction_id, 'astatine_mapped': astatine_smiles}
            writer.writerow(new_row)

def generate_mapped_data(input_file: str, output_file: str):
    """Generate a CSV file with atom-mapped reactions from astatine-mapped reactions."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['rxn_idx', 'atom_mapped']  # Assuming 'rxn_idx' is the ID column
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for row in reader:
            reaction_id = row['rxn_idx']  # Use the correct column name for reaction ID
            astatine_smiles = row['astatine_mapped']  # Use the correct column name for astatine SMILES
            try:
                atom_mapped_smiles = map_atoms(astatine_smiles)
                row = {'rxn_idx': reaction_id, 'atom_mapped': atom_mapped_smiles}
            except ValueError as e:
                print(f"Error processing reaction {reaction_id}: {e}")
                row = {'rxn_idx': reaction_id, 'atom_mapped': 'ERROR'}
            writer.writerow(row)

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    paths = load_paths('pipeline/config/paths.yaml')

    # Ensure all necessary directories exist
    ensure_directories(paths)

    # Step 1: Generate filtered reaction data
    generate_filtered_reaction_data(paths['raw_data'], paths['filtered_data'])
    print(f"Filtered reaction data generated and saved to {paths['filtered_data']}")

    # Step 2: Generate astatine data
    generate_astatine_data(paths['filtered_data'], paths['astatine_data'])
    print(f"Astatine test data generated and saved to {paths['astatine_data']}")

    # Step 3: Generate atom-mapped data
    generate_mapped_data(paths['astatine_data'], paths['mapped_data'])
    print(f"Atom-mapped test data generated and saved to {paths['mapped_data']}")

if __name__ == "__main__":
    main()
