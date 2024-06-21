import csv
import os
import sys
from evodex.mapping import map_atoms

# This deals with path issues
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def generate_mapped_data(input_file: str, output_file: str):
    """
    Generate a CSV file with atom-mapped reactions from astatine-mapped reactions.

    Parameters:
    input_file (str): Path to the input CSV file with astatine-mapped reactions.
    output_file (str): Path to the output CSV file to save atom-mapped reactions.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'atom_mapped']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in reader:
            reaction_id = row['id']
            astatine_smiles = row[reader.fieldnames[-1]]  # Last field
            try:
                atom_mapped_smiles = map_atoms(astatine_smiles)
                row = {'id': reaction_id, 'atom_mapped': atom_mapped_smiles}
            except ValueError as e:
                print(f"Error processing reaction {reaction_id}: {e}")
                row = {'id': reaction_id, 'atom_mapped': 'ERROR'}
            writer.writerow(row)

if __name__ == "__main__":
    input_file = 'tests/data/mapping_test_data.csv'
    output_file = 'tests/data/splitting_test_data.csv'
    generate_mapped_data(input_file, output_file)
    print(f"Atom-mapped test data generated and saved to {output_file}")
