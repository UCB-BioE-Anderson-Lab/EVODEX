import os
import pandas as pd
import json

def load_raw_data(data_dir):
    """Loads raw reactions and evodex data files from the specified directory."""
    print("Loading data...")
    raw_reactions = pd.read_csv(f'{data_dir}/raw_reactions.csv')
    evodex_files = {
        'C': 'EVODEX-C_reaction_operators.csv',
        'Cm': 'EVODEX-Cm_reaction_operators.csv',
        'E': 'EVODEX-E_reaction_operators.csv',
        'Em': 'EVODEX-Em_reaction_operators.csv',
        'N': 'EVODEX-N_reaction_operators.csv',
        'Nm': 'EVODEX-Nm_reaction_operators.csv',
        'R': 'EVODEX-R_full_reactions.csv',
        'P': 'EVODEX-P_partial_reactions.csv',
        'M': 'EVODEX-M_unique_masses.csv',
        'F': 'EVODEX-F_unique_formulas.csv'
    }

    evodex_data = {}
    for evodex_type, file_name in evodex_files.items():
        evodex_data[evodex_type] = pd.read_csv(f'{data_dir}/{file_name}')
        print(f'Loaded {evodex_type} data with {evodex_data[evodex_type].shape[0]} rows')
    
    return raw_reactions, evodex_data

def expand_sources(df, source_column='sources', delimiter=','):
    """Expands source column into separate rows for better processing."""
    df[source_column] = df[source_column].astype(str).str.replace('"', '')
    df_expanded = df.assign(**{source_column: df[source_column].str.split(delimiter)}).explode(source_column)
    return df_expanded

def map_evodex_to_ec(raw_reactions, evodex_data):
    ec_map = {'R': {}, 'P': {}, 'C': {}, 'Cm': {}, 'E': {}, 'Em': {}, 'N': {}, 'Nm': {}, 'M': {}, 'F': {}}
    
    # Ensure rxn_idx is treated as strings
    raw_reactions['rxn_idx'] = raw_reactions['rxn_idx'].astype(str)
    ec_dict = raw_reactions.set_index('rxn_idx')['ec_num'].to_dict()

    # Process R
    evodex_r_expanded = expand_sources(evodex_data['R'], 'sources')
    for _, row in evodex_r_expanded.iterrows():
        evodex_id = str(row['id'])
        source_ids = str(row['sources']).split(',')
        for source_id in source_ids:
            source_id = source_id.strip()
            ec_num = ec_dict.get(source_id, None)
            if ec_num:
                if evodex_id not in ec_map['R']:
                    ec_map['R'][evodex_id] = set()
                ec_map['R'][evodex_id].add(ec_num)
    print(f"Mapped EVODEX-R: {len(ec_map['R'])} entries")

    # Process P
    evodex_p_expanded = expand_sources(evodex_data['P'], 'sources')
    for _, row in evodex_p_expanded.iterrows():
        evodex_id = str(row['id'])
        source_ids = str(row['sources']).split(',')
        ec_set = set()
        for source_id in source_ids:
            source_id = source_id.strip()
            ec_set.update(ec_map['R'].get(source_id, set()))
        ec_map['P'][evodex_id] = ec_set
    print(f"Mapped EVODEX-P: {len(ec_map['P'])} entries")

    # Process other EVODEX types
    for evodex_type in ['C', 'Cm', 'E', 'Em', 'N', 'Nm', 'M', 'F']:
        evodex_expanded = expand_sources(evodex_data[evodex_type], 'sources')
        for _, row in evodex_expanded.iterrows():
            evodex_id = str(row['id'])
            source_ids = str(row['sources']).split(',')
            if evodex_id not in ec_map[evodex_type]:
                ec_map[evodex_type][evodex_id] = set()
            for source_id in source_ids:
                source_id = source_id.strip()
                ec_map[evodex_type][evodex_id].update(ec_map['P'].get(source_id, set()))
        print(f"Mapped EVODEX-{evodex_type}: {len(ec_map[evodex_type])} entries")

    return ec_map

def main():
    data_dir = 'website/data'
    output_dir = 'website/analysis'
    os.makedirs(output_dir, exist_ok=True)

    raw_reactions, evodex_data = load_raw_data(data_dir)
    ec_map = map_evodex_to_ec(raw_reactions, evodex_data)
    
    # Write ec_map to JSON file
    ec_map_serializable = {k: {kk: list(vv) for kk, vv in v.items()} for k, v in ec_map.items()}
    with open(f'{output_dir}/ec_map.json', 'w') as f:
        json.dump(ec_map_serializable, f, indent=4)
    print(f'EC map written to {output_dir}/ec_map.json')
    
    return ec_map

if __name__ == "__main__":
    main()
