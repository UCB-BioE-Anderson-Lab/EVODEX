import csv
import os
from evodex.decofactor import remove_cofactors
from evodex.operators import extract_operator
from evodex.formula import calculate_formula_diff, calculate_exact_mass
from pipeline.config import load_paths

def load_reactions(input_csv: str) -> list[dict]:
    reactions = []
    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        headers = reader.fieldnames
        print(f"Headers in {input_csv}: {headers}")  # Debugging statement
        for row in reader:
            reactions.append({'rxn_idx': row['rxn_idx'], 'atom_mapped': row['atom_mapped']})
    return reactions

def generate_partial_reactions(input_csv: str, output_csv: str):
    reactions = load_reactions(input_csv)
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['rxn_idx', 'partial_reaction_smiles'])

        for reaction in reactions:
            reaction_id = reaction['rxn_idx']
            reaction_smiles = reaction['atom_mapped']
            try:
                partial_reaction = remove_cofactors(reaction_smiles)
                writer.writerow([reaction_id, partial_reaction])
            except ValueError as e:
                print(f"No partial reactions generated for reaction ID {reaction_id}: {e}")
            except Exception as e:
                print(f"Error processing reaction ID {reaction_id}: {e}")

def load_partial_reactions(input_csv: str) -> list[dict]:
    reactions = []
    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        headers = reader.fieldnames
        print(f"Headers in {input_csv}: {headers}")  # Debugging statement
        for row in reader:
            reactions.append({'rxn_idx': row['rxn_idx'], 'partial_reaction_smiles': row.get('partial_reaction_smiles', row.get('operator_smirks'))})
    return reactions

def generate_operator_data(input_csv: str, output_csv: str, extract_params: dict):
    reactions = load_partial_reactions(input_csv)
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['rxn_idx', 'operator_smirks'])

        for reaction in reactions:
            reaction_id = reaction['rxn_idx']
            partial_reaction_smiles = reaction['partial_reaction_smiles']
            try:
                operator_smirks = extract_operator(partial_reaction_smiles, **extract_params)
                writer.writerow([reaction_id, operator_smirks])
            except ValueError as e:
                print(f"No operator SMIRKS generated for reaction ID {reaction_id}: {e}")
            except Exception as e:
                print(f"Error processing reaction ID {reaction_id}: {e}")

def generate_formula_mass_data(input_csv: str, output_csv: str):
    reactions = load_partial_reactions(input_csv)
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['rxn_idx', 'formula', 'mass'])

        for reaction in reactions:
            reaction_id = reaction['rxn_idx']
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

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    paths = load_paths('pipeline/config/paths.yaml')

    # Ensure all necessary directories exist
    ensure_directories(paths)

    # Step 1: Generate partial reactions
    generate_partial_reactions(paths['mapped_data'], paths['partial_reactions'])
    print(f"Partial reactions generated and saved to {paths['partial_reactions']}")
    
    # Step 2: Generate operator data
    extract_params = {
        'include_stereochemistry': True,
        'include_sigma': True,
        'include_pi': True,
        'include_unmapped_hydrogens': True,
        'include_unmapped_heavy_atoms': True,
        'include_static_hydrogens': False
    }
    generate_operator_data(paths['partial_reactions'], paths['operator_data'], extract_params)
    print(f"Operator data generated and saved to {paths['operator_data']}")
    
    # Step 3: Generate formula and mass data
    generate_formula_mass_data(paths['operator_data'], paths['final_mf_data'])
    print(f"Formula and mass data generated and saved to {paths['final_mf_data']}")

if __name__ == "__main__":
    main()
