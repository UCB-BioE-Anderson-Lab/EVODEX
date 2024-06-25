import csv
from evodex.operators import extract_operator

def load_partial_reactions(input_csv: str) -> list[dict]:
    reactions = []
    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append({'id': row['id'], 'partial_reaction_smiles': row['partial_reaction_smiles']})
    return reactions

def generate_operator_data(input_csv: str, output_csv: str, extract_params: dict):
    reactions = load_partial_reactions(input_csv)
    with open(output_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['id', 'operator_smirks'])

        for reaction in reactions:
            reaction_id = reaction['id']
            partial_reaction_smiles = reaction['partial_reaction_smiles']
            try:
                operator_smirks = extract_operator(partial_reaction_smiles, **extract_params)
                writer.writerow([reaction_id, operator_smirks])
            except ValueError as e:
                print(f"No operator SMIRKS generated for reaction ID {reaction_id}: {e}")
            except Exception as e:
                print(f"Error processing reaction ID {reaction_id}: {e}")

if __name__ == "__main__":
    input_csv = 'tests/data/operator_test_data.csv'
    output_csv = 'tests/data/final_ro_data.csv'

    # Define the parameters for extracting the operator SMIRKS
    extract_params = {
        'include_stereochemistry': True,
        'include_sigma': True,
        'include_pi': True,
        'include_unmapped_hydrogens': True,
        'include_unmapped_heavy_atoms': True,
        'include_static_hydrogens': False
    }
    
    # Generating each EVODEX file with different parameters
    ro_metadata = {
        'E': {
            'params': {
                'include_stereochemistry': True,
                'include_sigma': True,
                'include_pi': True,
                'include_unmapped_hydrogens': True,
                'include_unmapped_heavy_atoms': True,
                'include_static_hydrogens': False
            },
            'filename': 'EVODEX-E_reaction_operators.csv'
        },
        'C': {
            'params': {
                'include_stereochemistry': True,
                'include_sigma': False,
                'include_pi': False,
                'include_unmapped_hydrogens': True,
                'include_unmapped_heavy_atoms': True,
                'include_static_hydrogens': False
            },
            'filename': 'EVODEX-C_reaction_operators.csv'
        },
        'N': {
            'params': {
                'include_stereochemistry': True,
                'include_sigma': True,
                'include_pi': False,
                'include_unmapped_hydrogens': True,
                'include_unmapped_heavy_atoms': True,
                'include_static_hydrogens': False
            },
            'filename': 'EVODEX-N_reaction_operators.csv'
        },
        'Em': {
            'params': {
                'include_stereochemistry': False,
                'include_sigma': True,
                'include_pi': True,
                'include_unmapped_hydrogens': False,
                'include_unmapped_heavy_atoms': False,
                'include_static_hydrogens': False
            },
            'filename': 'EVODEX-Em_reaction_operators.csv'
        },
        'Cm': {
            'params': {
                'include_stereochemistry': False,
                'include_sigma': False,
                'include_pi': False,
                'include_unmapped_hydrogens': False,
                'include_unmapped_heavy_atoms': False,
                'include_static_hydrogens': False
            },
            'filename': 'EVODEX-Cm_reaction_operators.csv'
        },
        'Nm': {
            'params': {
                'include_stereochemistry': False,
                'include_sigma': True,
                'include_pi': False,
                'include_unmapped_hydrogens': False,
                'include_unmapped_heavy_atoms': False,
                'include_static_hydrogens': False
            },
            'filename': 'EVODEX-Nm_reaction_operators.csv'
        }
    }

    for evodex_type, metadata in ro_metadata.items():
        output_csv = metadata['filename']
        extract_params = metadata['params']
        generate_operator_data(input_csv, output_csv, extract_params)
        print(f"Operator data generated and saved to {output_csv}")
