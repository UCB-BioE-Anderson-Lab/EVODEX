import csv
from evodex.splitting import split_reaction

def load_reactions(input_csv: str) -> list[dict]:
    reactions = []
    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append(row)
    return reactions

def pick_smallest_partial(reaction_smiles: list[str]) -> str:
    smallest_partial = min(reaction_smiles, key=lambda x: len(x.split('>>')[0].split('.')))
    return smallest_partial

def generate_partial_reactions(input_csv: str, output_csv: str):
    reactions = load_reactions(input_csv)
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['id', 'partial_reaction_smiles'])

        for reaction in reactions:
            reaction_id = reaction['id']
            reaction_smiles = reaction['atom_mapped']
            partial_reactions = split_reaction(reaction_smiles)
            if partial_reactions:
                smallest_partial = pick_smallest_partial(partial_reactions)
                writer.writerow([reaction_id, smallest_partial])
            else:
                print(f"No partial reactions generated for reaction ID {reaction_id}")

if __name__ == "__main__":
    input_csv = 'tests/data/splitting_test_data.csv'
    output_csv = 'tests/data/operator_test_data.csv'
    generate_partial_reactions(input_csv, output_csv)
    print(f"Partial reactions generated and saved to {output_csv}")
