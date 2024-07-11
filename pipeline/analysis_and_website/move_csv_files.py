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
    error_logs = []

    for filename in os.listdir(source_dir):
        if filename.endswith(".csv"):
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

    raw_src_path = paths['raw_data']
    raw_dest_path = os.path.join(data_dir, os.path.basename(raw_src_path))
    
    # Copy raw data file without conversion
    shutil.copy(raw_src_path, raw_dest_path)
    print(f"Copied: {raw_src_path} to {raw_dest_path}")

    # Save error logs
    if error_logs:
        error_log_df = pd.DataFrame(error_logs)
        error_log_df.to_csv(os.path.join(error_dir, 'astatine_conversion_errors.csv'), index=False)
        print(f"Conversion errors logged to {os.path.join(data_dir, 'astatine_conversion_errors.csv')}")

if __name__ == "__main__":
    data_dir = 'website/data'
    paths = load_paths('pipeline/config/paths.yaml')
    move_and_convert_csv_files(data_dir, paths)
