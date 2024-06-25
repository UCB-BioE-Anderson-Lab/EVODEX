import csv
import os
import sys
from evodex.filter import validate_smiles

# This deals with path issues
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def generate_filtered_reaction_data(input_file: str, output_file: str):
    """
    Generate a CSV file with validated reactions from the input file.

    Parameters:
    input_file (str): Path to the input CSV file with reactions.
    output_file (str): Path to the output CSV file to save validated reactions.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in reader:
            smiles = row['rxn']
            if validate_smiles(smiles):
                writer.writerow(row)

if __name__ == "__main__":
    input_file = 'tests/data/raw_test_reactions.csv'
    output_file = 'tests/data/astatine_test_data.csv'
    generate_filtered_reaction_data(input_file, output_file)
    print(f"Filtered reaction data generated and saved to {output_file}")
