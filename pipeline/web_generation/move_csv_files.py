import os
import shutil
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pipeline.config import load_paths
from evodex.astatine import astatine_to_hydrogen_reaction

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

def move_and_convert_csv_files(data_dir, paths):
    source_dir = os.path.join('data', 'processed')
    error_dir = os.path.join('data', 'errors')
    evodex_dir = 'evodex/data'
    error_logs = []

    # Move and convert other CSV files
    for filename in os.listdir(source_dir):
        if filename.endswith(".csv") and filename != 'EVODEX-E_synthesis_subset.csv':
            src_path = os.path.join(source_dir, filename)
            dest_path = os.path.join(data_dir, filename)

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
    synthesis_dest_path = os.path.join(data_dir, 'EVODEX-E_synthesis_subset.csv')

    if os.path.exists(synthesis_src_path):
        synthesis_df = pd.read_csv(synthesis_src_path)

        # Read the already converted and copied EVODEX-E file
        evodex_e_dest_path = os.path.join(data_dir, os.path.basename(paths['evodex_e']))
        evodex_e_df = pd.read_csv(evodex_e_dest_path)

        # Create a map for EVODEX-E ID to SMIRKS
        evodex_e_smirks_map = evodex_e_df.set_index('id')['smirks'].to_dict()

        # Add SMIRKS column to the synthesis subset DataFrame
        synthesis_df['smirks'] = synthesis_df['id'].map(evodex_e_smirks_map)

        # Save the modified DataFrame to the destination path
        synthesis_df.to_csv(synthesis_dest_path, index=False)
        print(f"Copied and converted: {synthesis_src_path} to {synthesis_dest_path}")

    selected_src_path = paths['selected_data']
    selected_dest_path = os.path.join(data_dir, 'selected_data.csv')
    
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
        src_path = os.path.join(data_dir, filename)
        dest_path = os.path.join(evodex_dir, filename)
        if os.path.exists(src_path):
            shutil.copy(src_path, dest_path)
            print(f"Copied: {src_path} to {dest_path}")

if __name__ == "__main__":
    data_dir = 'website/data'
    paths = load_paths('pipeline/config/paths.yaml')
    move_and_convert_csv_files(data_dir, paths)
