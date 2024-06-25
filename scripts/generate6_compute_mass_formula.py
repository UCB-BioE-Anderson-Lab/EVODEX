import csv
from evodex.operators import extract_operator
from evodex.formula import calculate_formula_diff, calculate_exact_mass

def load_partial_reactions(input_csv: str) -> list[dict]:
    reactions = []
    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append({'id': row['id'], 'partial_reaction_smiles': row['partial_reaction_smiles']})
    return reactions

def generate_formula_mass_data(input_csv: str, output_csv: str):
    reactions = load_partial_reactions(input_csv)
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['rxn_id', 'formula', 'mass'])

        for reaction in reactions:
            reaction_id = reaction['id']
            partial_reaction_smiles = reaction['partial_reaction_smiles']
            try:
                # Calculate the formula difference
                formula_diff = calculate_formula_diff(partial_reaction_smiles)
                # Convert the formula_diff to a string key for mapping
                formula_key = str(sorted(formula_diff.items()))
                # Calculate the exact mass
                mass_difference = calculate_exact_mass(formula_diff)
                writer.writerow([reaction_id, formula_key, mass_difference])
            except ValueError as e:
                print(f"Error calculating formula or mass for reaction ID {reaction_id}: {e}")
            except Exception as e:
                print(f"Error processing reaction ID {reaction_id}: {e}")

if __name__ == "__main__":
    input_csv = 'tests/data/operator_test_data.csv'
    output_csv = 'tests/data/final_mf_data.csv'
    generate_formula_mass_data(input_csv, output_csv)
    print(f"Formula and mass data generated and saved to {output_csv}")
