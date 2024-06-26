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

def generate_full_reactions(input_csv: str, output_csv: str):
    reactions = load_reactions(input_csv)
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['rxn_idx', 'full_reaction'])

        for reaction in reactions:
            reaction_id = reaction['rxn_idx']
            reaction_smiles = reaction['atom_mapped']
            try:
                full_reaction = extract_operator(reaction_smiles, include_stereochemistry=True)
                writer.writerow([reaction_id, full_reaction])
            except ValueError as e:
                print(f"No full reaction generated for reaction ID {reaction_id}: {e}")
            except Exception as e:
                print(f"Error processing reaction ID {reaction_id}: {e}")

def generate_partial_reactions(input_csv: str, output_csv: str):
    reactions = []
    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append({'rxn_idx': row['rxn_idx'], 'full_reaction': row['full_reaction']})
    
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['rxn_idx', 'partial_reaction_smiles'])

        for reaction in reactions:
            reaction_id = reaction['rxn_idx']
            reaction_smiles = reaction['full_reaction']
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
        writer.writerow(['rxn_idx', 'operator_smirks', 'id'])

        for reaction in reactions:
            reaction_id = reaction['rxn_idx']
            partial_reaction_smiles = reaction['partial_reaction_smiles']
            try:
                operator_smirks = extract_operator(partial_reaction_smiles, **extract_params)
                writer.writerow([reaction_id, operator_smirks, reaction_id])
            except ValueError as e:
                print(f"No operator SMIRKS generated for reaction ID {reaction_id}: {e}")
            except Exception as e:
                print(f"Error processing reaction ID {reaction_id}: {e}")

def generate_all_operator_data(input_csv: str, paths: dict):
    extract_params_map = {
        'evodex_e': {
            'include_stereochemistry': True,
            'include_sigma': True,
            'include_pi': True,
            'include_unmapped_hydrogens': True,
            'include_unmapped_heavy_atoms': True,
            'include_static_hydrogens': False
        },
        'evodex_c': {
            'include_stereochemistry': True,
            'include_sigma': False,
            'include_pi': False,
            'include_unmapped_hydrogens': True,
            'include_unmapped_heavy_atoms': True,
            'include_static_hydrogens': False
        },
        'evodex_n': {
            'include_stereochemistry': True,
            'include_sigma': True,
            'include_pi': False,
            'include_unmapped_hydrogens': True,
            'include_unmapped_heavy_atoms': True,
            'include_static_hydrogens': False
        },
        'evodex_em': {
            'include_stereochemistry': False,
            'include_sigma': True,
            'include_pi': True,
            'include_unmapped_hydrogens': False,
            'include_unmapped_heavy_atoms': False,
            'include_static_hydrogens': False
        },
        'evodex_cm': {
            'include_stereochemistry': False,
            'include_sigma': False,
            'include_pi': False,
            'include_unmapped_hydrogens': False,
            'include_unmapped_heavy_atoms': False,
            'include_static_hydrogens': False
        },
        'evodex_nm': {
            'include_stereochemistry': False,
            'include_sigma': True,
            'include_pi': False,
            'include_unmapped_hydrogens': False,
            'include_unmapped_heavy_atoms': False,
            'include_static_hydrogens': False
        }
    }

    for key, params in extract_params_map.items():
        generate_operator_data(input_csv, paths[key], params)
        print(f"Operator data generated and saved to {paths[key]}")

def generate_formula_mass_data(input_csv: str, output_csv_f: str, output_csv_m: str):
    reactions = load_partial_reactions(input_csv)
    with open(output_csv_f, 'w', newline='') as file_f, open(output_csv_m, 'w', newline='') as file_m:
        writer_f = csv.writer(file_f)
        writer_m = csv.writer(file_m)
        writer_f.writerow(['rxn_idx', 'formula'])
        writer_m.writerow(['rxn_idx', 'mass'])

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
                writer_f.writerow([reaction_id, formula_key])
                writer_m.writerow([reaction_id, mass_difference])
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

    # Step 1: Generate full reactions
    generate_full_reactions(paths['mapped_data'], paths['evodex_r'])
    print(f"Full reactions generated and saved to {paths['evodex_r']}")

    # Step 2: Generate partial reactions
    generate_partial_reactions(paths['evodex_r'], paths['evodex_p'])
    print(f"Partial reactions generated and saved to {paths['evodex_p']}")
    
    # Step 3: Generate all operator data
    generate_all_operator_data(paths['evodex_p'], paths)
    
    # Step 4: Generate formula and mass data
    generate_formula_mass_data(paths['evodex_e'], paths['evodex_f'], paths['evodex_m'])
    print(f"Formula and mass data generated and saved to {paths['evodex_f']} and {paths['evodex_m']}")

if __name__ == "__main__":
    main()
