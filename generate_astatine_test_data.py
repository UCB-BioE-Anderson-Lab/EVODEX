import csv
import os
import sys
from evodex.astatine import hydrogen_to_astatine

# This deals with path issues
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def generate_astatine_data(input_file: str, output_file: str):
    """
    Generate a CSV file with astatine-containing reactions from hydrogen-containing reactions.

    Parameters:
    input_file (str): Path to the input CSV file with hydrogen-containing reactions.
    output_file (str): Path to the output CSV file to save astatine-containing reactions.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['astatine_mapped']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in reader:
            hydrogen_smiles = row['rxn']
            astatine_smiles = hydrogen_to_astatine(hydrogen_smiles)
            row['astatine_mapped'] = astatine_smiles
            writer.writerow(row)

if __name__ == "__main__":
    input_file = 'tests/data/test_reactions_with_hydrogen.csv'
    output_file = 'tests/data/test_reactions_with_astatine.csv'
    generate_astatine_data(input_file, output_file)
    print(f"Astatine test data generated and saved to {output_file}")