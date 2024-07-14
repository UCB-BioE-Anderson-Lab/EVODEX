import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive plots
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency

def load_data(data_dir):
    """
    Load raw reactions and EVODEX data files.
    """
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
        'P': 'EVODEX-P_partial_reactions.csv'
    }

    evodex_data = {}
    for evodex_type, file_name in evodex_files.items():
        evodex_data[evodex_type] = pd.read_csv(f'{data_dir}/{file_name}')
        print(f'Loaded {evodex_type} data with {evodex_data[evodex_type].shape[0]} rows')
    
    return raw_reactions, evodex_data

def expand_sources(df, source_column='sources', delimiter=','):
    """
    Expand the sources column into separate rows.
    """
    df[source_column] = df[source_column].astype(str).str.replace('"', '')
    df_expanded = df.assign(**{source_column: df[source_column].str.split(delimiter)}).explode(source_column)
    return df_expanded

def split_ec_number(ec_num):
    """
    Split EC number into its component levels.
    """
    parts = ec_num.split('.')
    return parts + [''] * (4 - len(parts))

def map_ec_numbers(raw_reactions, evodex_data):
    """
    Map EC numbers to EVODEX types.
    """
    ec_map = raw_reactions[['rxn_idx', 'ec_num']].dropna().copy()
    ec_map['rxn_idx'] = ec_map['rxn_idx'].astype(str)
    print("EC Map (rxn_idx to ec_num):")
    print(ec_map.head())
    
    # Expand EVODEX-R to raw reactions
    evodex_r = expand_sources(evodex_data['R'], source_column='sources', delimiter=',')
    evodex_r = evodex_r.merge(ec_map, left_on='sources', right_on='rxn_idx', how='left')
    print("EVODEX-R with EC numbers:")
    print(evodex_r.head())
    
    # Expand EVODEX-P to EVODEX-R
    evodex_p = expand_sources(evodex_data['P'], source_column='sources', delimiter=',')
    evodex_p = evodex_p.merge(evodex_r[['id', 'ec_num']], left_on='sources', right_on='id', how='left', suffixes=('', '_r'))
    print("EVODEX-P with EC numbers:")
    print(evodex_p.head())
    
    # Map other EVODEX types to EVODEX-P
    mapped_data = pd.DataFrame()
    for evodex_type in ['C', 'Cm', 'E', 'Em', 'N', 'Nm']:
        print(f'Processing {evodex_type}...')
        evodex_df = expand_sources(evodex_data[evodex_type], source_column='sources', delimiter=',')
        evodex_df = evodex_df.merge(evodex_p[['id', 'ec_num']], left_on='sources', right_on='id', how='left', suffixes=('', '_p'))
        evodex_df['evodex_type'] = evodex_type
        mapped_data = pd.concat([mapped_data, evodex_df], axis=0)
        print(f"{evodex_type} mapping complete. Data snapshot:")
        print(evodex_df.head())
    
    print("Mapped data (first few rows):")
    print(mapped_data.head())
    return mapped_data

def analyze_ec_levels(mapped_data, output_dir):
    """
    Analyze and create heatmaps for EC levels.
    """
    print("Analyzing EC levels...")
    mapped_data[['ec_level_1', 'ec_level_2', 'ec_level_3', 'ec_level_4']] = mapped_data['ec_num'].apply(lambda x: pd.Series(split_ec_number(x)))
    
    for level in range(1, 5):
        ec_level = f'ec_level_{level}'
        ec_distribution = mapped_data.groupby(['evodex_type', ec_level]).size().reset_index(name='count')
        print(f'EC Level {level} distribution:')
        print(ec_distribution)
        
        ec_pivot = ec_distribution.pivot(index=ec_level, columns='evodex_type', values='count').fillna(0)
        
        plt.figure(figsize=(12, 6))
        sns.heatmap(ec_pivot, cmap='coolwarm', annot=True, fmt='g')
        plt.title(f'Correlation of EC Level {level} Numbers with EVODEX Types')
        plt.savefig(os.path.join(output_dir, f'ec_level_{level}_correlation_heatmap.png'))
        ec_distribution.to_csv(os.path.join(output_dir, f'ec_level_{level}_distribution.csv'), index=False)
        
        # Chi-square test
        chi2, p, dof, ex = chi2_contingency(ec_pivot)
        with open(os.path.join(output_dir, f'ec_level_{level}_correlation_stats.txt'), 'w') as f:
            f.write(f'Chi-square test\nChi2: {chi2}\nP-value: {p}\nDegrees of Freedom: {dof}\nExpected Frequencies: {ex}\n')
        print(f'Creating heatmap for EC level {level}')

def main():
    data_dir = 'website/data'
    output_dir = 'website/pages'
    os.makedirs(output_dir, exist_ok=True)
    
    raw_reactions, evodex_data = load_data(data_dir)
    mapped_data = map_ec_numbers(raw_reactions, evodex_data)
    analyze_ec_levels(mapped_data, output_dir)
    print("EC Number Correlation analysis complete.")

if __name__ == "__main__":
    main()
