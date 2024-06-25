import csv
from evodex.decofactor import remove_cofactors

def load_reactions(input_csv: str) -> list[dict]:
    reactions = []
    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append({'id': row['id'], 'atom_mapped': row['atom_mapped']})
    return reactions

def generate_partial_reactions(input_csv: str, output_csv: str):
    reactions = load_reactions(input_csv)
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['id', 'partial_reaction_smiles'])

        for reaction in reactions:
            reaction_id = reaction['id']
            reaction_smiles = reaction['atom_mapped']
            try:
                partial_reaction = remove_cofactors(reaction_smiles)
                writer.writerow([reaction_id, partial_reaction])
            except ValueError as e:
                print(f"No partial reactions generated for reaction ID {reaction_id}: {e}")
            except Exception as e:
                print(f"Error processing reaction ID {reaction_id}: {e}")

if __name__ == "__main__":
    input_csv = 'tests/data/splitting_test_data.csv'
    output_csv = 'tests/data/operator_test_data.csv'
    generate_partial_reactions(input_csv, output_csv)
    print(f"Partial reactions generated and saved to {output_csv}")
