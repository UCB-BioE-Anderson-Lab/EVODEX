import os
import shutil
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pipeline.config import load_paths
from evodex.astatine import astatine_to_hydrogen_reaction
from ast import literal_eval
import json

def convert_smiles_column(df, column_name):
    errors = []
    def convert_smiles(smiles):
        try:
            return astatine_to_hydrogen_reaction(smiles) if pd.notnull(smiles) else smiles
        except Exception as e:
            errors.append((smiles, str(e)))
            return smiles
    
    df[column_name] = df[column_name].apply(convert_smiles)
    return df, errors

def move_and_convert_csv_files(evodex_dir, paths):
    source_dir = os.path.join('data', 'processed')
    error_dir = os.path.join('data', 'errors')
    error_logs = []

    # Move and convert other CSV files
    for filename in os.listdir(source_dir):
        if filename.endswith(".csv") and filename != 'EVODEX-E_synthesis_subset.csv':
            src_path = os.path.join(source_dir, filename)
            dest_path = os.path.join(evodex_dir, filename)

            # Read the CSV file into a DataFrame
            df = pd.read_csv(src_path)

            # Convert the SMILES/SMIRKS column if it exists
            if 'smirks' in df.columns:
                df, errors = convert_smiles_column(df, 'smirks')
                for smiles, error in errors:
                    error_logs.append({'id': df[df['smirks'] == smiles]['id'].values[0], 'smirks': smiles, 'sources': df[df['smirks'] == smiles]['sources'].values[0], 'error_message': error})
            if 'smiles' in df.columns:
                df, errors = convert_smiles_column(df, 'smiles')
                for smiles, error in errors:
                    error_logs.append({'id': df[df['smiles'] == smiles]['id'].values[0], 'smiles': smiles, 'sources': df[df['smiles'] == smiles]['sources'].values[0], 'error_message': error})

            # Save the modified DataFrame to the destination path
            df.to_csv(dest_path, index=False)
            print(f"Copied and converted: {src_path} to {dest_path}")

    # Handle EVODEX-E_synthesis_subset.csv
    synthesis_src_path = os.path.join(source_dir, 'EVODEX-E_synthesis_subset.csv')
    synthesis_dest_path = os.path.join(evodex_dir, 'EVODEX-E_synthesis_subset.csv')

    if os.path.exists(synthesis_src_path):
        synthesis_df = pd.read_csv(synthesis_src_path)

        # Read the already converted and copied EVODEX-E file
        evodex_e_dest_path = os.path.join(evodex_dir, os.path.basename(paths['evodex_e']))
        evodex_e_df = pd.read_csv(evodex_e_dest_path)

        # Create a map for EVODEX-E ID to SMIRKS
        evodex_e_smirks_map = evodex_e_df.set_index('id')['smirks'].to_dict()

        # Add SMIRKS column to the synthesis subset DataFrame
        synthesis_df['smirks'] = synthesis_df['id'].map(evodex_e_smirks_map)

        # Save the modified DataFrame to the destination path
        synthesis_df.to_csv(synthesis_dest_path, index=False)
        print(f"Copied and converted: {synthesis_src_path} to {synthesis_dest_path}")

    selected_src_path = paths['selected_data']
    selected_dest_path = os.path.join(evodex_dir, 'selected_data.csv')
    
    # Copy selected data file without conversion
    shutil.copy(selected_src_path, selected_dest_path)
    print(f"Copied: {selected_src_path} to {selected_dest_path}")

    # Save error logs
    if error_logs:
        os.makedirs(error_dir, exist_ok=True)
        error_log_df = pd.DataFrame(error_logs)
        error_log_df.to_csv(os.path.join(error_dir, 'astatine_conversion_errors.csv'), index=False)
        print(f"Conversion errors logged to {os.path.join(error_dir, 'astatine_conversion_errors.csv')}")

    # Copy specific files from EVODEX/website to EVODEX/evodex/data
    files_to_copy = [
        'EVODEX-C_reaction_operators.csv', 'EVODEX-N_reaction_operators.csv',
        'EVODEX-Cm_reaction_operators.csv', 'EVODEX-Nm_reaction_operators.csv',
        'EVODEX-E_reaction_operators.csv', 'EVODEX-E_synthesis_subset.csv',
        'EVODEX-Em_reaction_operators.csv', 'EVODEX-F_unique_formulas.csv',
        'EVODEX-M_mass_spec_subset.csv', 'EVODEX-M_unique_masses.csv'
    ]

    for filename in files_to_copy:
        src_path = os.path.join(evodex_dir, filename)
        dest_path = os.path.join(evodex_dir, filename)
        if os.path.exists(src_path) and src_path != dest_path:
            shutil.copy(src_path, dest_path)
            print(f"Copied: {src_path} to {dest_path}")

    # Generate optimized JSON structure
    generate_optimized_json(evodex_dir, os.path.join('website', 'data', 'optimized_data.json'))

def _parse_sources(sources):
    if pd.isnull(sources):
        return []
    if isinstance(sources, int):
        return [str(sources)]
    if isinstance(sources, str):
        return sources.replace('"', '').split(',')
    if isinstance(sources, list):
        return [str(source) for source in sources]
    return []

def generate_optimized_json(evodex_dir, output_json_path):
    data = {
        "EVODEX-R": {},
        "EVODEX-P": {},
        "EVODEX-F": {},
        "EVODEX-M": {},
        "EVODEX-C": {},
        "EVODEX-Cm": {},
        "EVODEX-E": {},
        "EVODEX-Em": {},
        "EVODEX-N": {},
        "EVODEX-Nm": {},
        "synthesis_subset": {},
        "mass_spec_subset": {}
    }

    def parse_sources(sources):
        if pd.isnull(sources):
            return []
        if isinstance(sources, int):
            return [str(sources)]
        if isinstance(sources, str):
            return sources.replace('"', '').split(',')
        if isinstance(sources, list):
            return [str(source) for source in sources]
        return []

    # Load selected_reactions.csv
    selected_reactions_path = os.path.join(evodex_dir, 'selected_reactions.csv')
    selected_reactions_df = pd.read_csv(selected_reactions_path)
    selected_reactions_map = selected_reactions_df.set_index('rxn_idx').to_dict(orient='index')

    # Create a dictionary to map EVODEX-P IDs to their children
    evodex_p_children_map = {}

    # Read and process CSV files to populate the JSON structure
    for filename in os.listdir(evodex_dir):
        if not filename.endswith(".csv") or not filename.startswith("EVODEX-"):
            print('Skipping:', filename)
            continue

        file_path = os.path.join(evodex_dir, filename)
        df = pd.read_csv(file_path)

        if 'smirks' in df.columns and 'sources' in df.columns and 'id' in df.columns:
            for _, row in df.iterrows():
                reaction_id = row['id']
                sources = parse_sources(row['sources'])
                smirks = row['smirks']

                # Handle specific EVODEX types
                if filename.startswith('EVODEX-R'):
                    details = {}
                    for source in sources:
                        source_int = int(source)
                        if source_int in selected_reactions_map:
                            details.update(selected_reactions_map[source_int])
                    
                    data["EVODEX-R"][reaction_id] = {
                        "id": reaction_id,
                        "smirks": smirks,
                        "sources": sources,
                        "details": {
                            "natural": details.get('natural', False),
                            "organism": details.get('organism', ''),
                            "protein_refs": details.get('protein_refs', '').split(',') if 'protein_refs' in details and isinstance(details.get('protein_refs', ''), str) and pd.notnull(details.get('protein_refs', '')) else [],
                            "protein_db": details.get('protein_db', ''),
                            "ec_num": details.get('ec_num', '')
                        }
                    }
                elif filename.startswith('EVODEX-P'):
                    data["EVODEX-P"][reaction_id] = {
                        "id": reaction_id,
                        "smirks": smirks,
                        "sources": sources,
                        "children": {}  # Initialize children dictionary
                    }
                elif filename.startswith(('EVODEX-E', 'EVODEX-N', 'EVODEX-C', 'EVODEX-Em', 'EVODEX-Nm', 'EVODEX-Cm')):
                    evodex_type = filename.split('_')[0]
                    if evodex_type not in data:
                        data[evodex_type] = {}
                    data[evodex_type][reaction_id] = {
                        "id": reaction_id,
                        "smirks": smirks,
                        "sources": sources
                    }

                    # Map the EVODEX-P to its children
                    for source in sources:
                        if source not in evodex_p_children_map:
                            evodex_p_children_map[source] = {}
                        if evodex_type not in evodex_p_children_map[source]:
                            evodex_p_children_map[source][evodex_type] = []
                        evodex_p_children_map[source][evodex_type].append(reaction_id)

        elif filename.startswith('EVODEX-F'):
            for _, row in df.iterrows():
                reaction_id = row['id']
                sources = parse_sources(row['sources'])
                formula = literal_eval(row['formula']) if 'formula' in row and isinstance(row['formula'], str) and pd.notnull(row['formula']) else ''
                
                data["EVODEX-F"][reaction_id] = {
                    "id": reaction_id,
                    "formula": formula,
                    "sources": sources
                }

                # Map the EVODEX-P to its children
                for source in sources:
                    if source not in evodex_p_children_map:
                        evodex_p_children_map[source] = {}
                    if 'EVODEX-F' not in evodex_p_children_map[source]:
                        evodex_p_children_map[source]['EVODEX-F'] = []
                    evodex_p_children_map[source]['EVODEX-F'].append(reaction_id)

        elif filename.startswith('EVODEX-M'):
            for _, row in df.iterrows():
                reaction_id = row['id']
                sources = parse_sources(row['sources'])
                mass = row['mass'] if 'mass' in row and pd.notnull(row['mass']) else ''

                data["EVODEX-M"][reaction_id] = {
                    "id": reaction_id,
                    "mass": mass,
                    "sources": sources
                }

                # Map the EVODEX-P to its children
                for source in sources:
                    if source not in evodex_p_children_map:
                        evodex_p_children_map[source] = {}
                    if 'EVODEX-M' not in evodex_p_children_map[source]:
                        evodex_p_children_map[source]['EVODEX-M'] = []
                    evodex_p_children_map[source]['EVODEX-M'].append(reaction_id)

        else: 
            raise ValueError('This should not happen:', filename)

    # Update EVODEX-P entries to include their children
    for p_id, children in evodex_p_children_map.items():
        if p_id in data["EVODEX-P"]:
            data["EVODEX-P"][p_id]["children"] = children

    # Save the JSON structure to a file
    with open(output_json_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)
    print(f"Optimized JSON saved to {output_json_path}")

if __name__ == "__main__":
    evodex_dir = 'evodex/data'
    paths = load_paths('pipeline/config/paths.yaml')
    move_and_convert_csv_files(evodex_dir, paths)
