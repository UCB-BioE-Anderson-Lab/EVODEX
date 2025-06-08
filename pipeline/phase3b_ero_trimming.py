import shutil
import os
import csv
import time
import pandas as pd
from pipeline.config import load_paths
from evodex.evaluation import match_operators
from evodex.astatine import convert_dataframe_smiles_column
from pipeline.version import __version__
import sys
csv.field_size_limit(sys.maxsize)

"""
Phase 3b: ERO Trimming (Dominance Pruning)

Goal: Perform dominance pruning of EVODEX-E, and trim P, F, and R accordingly.

Input:
- evodex_e_phase3a_preliminary (At)
- evodex_p_phase3a_retained (At)
- evodex_f_filtered (no SMIRKS, no conversion needed)
- evodex_r (At)

Outputs (all in data/processed/):
- evodex_e_phase3b_final (At)
- evodex_p_phase3b_final (At)
- evodex_f_phase3b_final (no SMIRKS, no conversion)
- evodex_r_phase3b_final (At)
"""

def main():
    start_time = time.time()
    print("Phase 3b ERO trimming (dominance pruning) started...")
    # Load paths
    paths = load_paths('pipeline/config/paths.yaml')
    
    # === Step 1: Convert required inputs to H ===
    print("Converting EVODEX-E Phase 3a preliminary to H...")
    evodex_e_df = pd.read_csv(paths['evodex_e_phase3a_preliminary'])
    converted_e_df, conversion_errors_e = convert_dataframe_smiles_column(evodex_e_df, 'smirks')
    converted_e_df['original_hash'] = evodex_e_df['id']
    converted_e_df['id'] = [f"EVODEX.{__version__}-E{idx+1}_temp" for idx in range(len(converted_e_df))]
    converted_e_df.to_csv(paths['evodex_e_phase3a_preliminary_H'], index=False)
    print(f"Saved: {paths['evodex_e_phase3a_preliminary_H']}")
    
    print("Converting EVODEX-P Phase 3a retained to H...")
    evodex_p_df = pd.read_csv(paths['evodex_p_phase3a_retained'])
    converted_p_df, conversion_errors_p = convert_dataframe_smiles_column(evodex_p_df, 'smirks')
    converted_p_df.to_csv(paths['evodex_p_phase3a_retained_H'], index=False)
    print(f"Saved: {paths['evodex_p_phase3a_retained_H']}")

    # === Step 2: Copy H versions to evodex/data for match_operators ===
    dst_e = os.path.join('evodex', 'data', 'EVODEX-E_reaction_operators.csv')
    shutil.copyfile(paths['evodex_e_phase3a_preliminary_H'], dst_e)
    print(f"Copied to {dst_e}")


    # Copy EVODEX-F to evodex/data
    dst_f = os.path.join('evodex', 'data', 'EVODEX-F_unique_formulas.csv')
    shutil.copyfile(paths['evodex_f_filtered'], dst_f)
    print(f"Copied to {dst_f}")

    # === Step 3: Run match_operators to build dominance table ===
    print("Running match_operators to build dominance table...")
    evodex_p_h_df = pd.read_csv(paths['evodex_p_phase3a_retained_H'])

    match_table = dict()
    for _, row in evodex_p_h_df.iterrows():
        p_id = row['id']
        smirks = row['smirks']
        matched_ops = match_operators(smirks)
        for op_id in matched_ops:
            if op_id not in match_table:
                match_table[op_id] = set()
            match_table[op_id].add(p_id)

    print(f"Total operators matched: {len(match_table)}")
    total_matches = sum(len(pids) for pids in match_table.values())
    print(f"Total matches recorded: {total_matches}")

    # === Extra statistics ===
    total_e_in_file = len(evodex_e_df)
    matched_e_ids = set(match_table.keys())
    unmatched_e_count = total_e_in_file - len(matched_e_ids)

    print(f"Total EVODEX-E operators in file: {total_e_in_file}")
    print(f"Number of operators with matches: {len(matched_e_ids)}")
    print(f"Number of operators with zero matches: {unmatched_e_count}")

    # Sample unmatched E IDs
    unmatched_e_df = evodex_e_df[~evodex_e_df['id'].isin(matched_e_ids)]
    if not unmatched_e_df.empty:
        print("Example unmatched operators (up to 5):")
        for _, row in unmatched_e_df.head(5).iterrows():
            print(f"ID: {row['id']}, SMIRKS: {row['smirks']}, Sources: {row['sources']}")

    # === Step 4: Perform dominance pruning ===
    operator_match_counts = {op_id: len(p_ids) for op_id, p_ids in match_table.items()}
    sorted_ops = sorted(operator_match_counts.items(), key=lambda x: x[1], reverse=True)

    retained_operators = set()
    explained_p = set()
    for op_id, count in sorted_ops:
        p_ids = match_table[op_id]
        if not p_ids.issubset(explained_p):
            retained_operators.add(op_id)
            explained_p.update(p_ids)

    print(f"Retained operators after dominance pruning: {len(retained_operators)}")

    total_p_in_file = len(evodex_p_h_df)
    matched_p_ids = explained_p
    unmatched_p_count = total_p_in_file - len(matched_p_ids)

    print(f"Total EVODEX-P IDs in file: {total_p_in_file}")
    print(f"Number of EVODEX-P IDs explained by at least one operator: {len(matched_p_ids)}")
    print(f"Number of unexplained (orphan) EVODEX-P IDs: {unmatched_p_count}")

    unmatched_p_df = evodex_p_h_df[~evodex_p_h_df['id'].isin(matched_p_ids)]
    if not unmatched_p_df.empty:
        print("Example unexplained EVODEX-P IDs (up to 5):")
        for _, row in unmatched_p_df.head(5).iterrows():
            print(row['id'])

    # Step 5: Gather surviving P hashes based on sources of retained E operators
    # Correct: Go back to original E dataframe to extract P hashes
    evodex_data_e_df = pd.read_csv('evodex/data/EVODEX-E_reaction_operators.csv')
    temp_to_hash_map = dict(zip(evodex_data_e_df['id'], evodex_data_e_df['original_hash']))
    retained_hashes = {temp_to_hash_map[op_id] for op_id in retained_operators}
    original_e_df = pd.read_csv(paths['evodex_e_phase3a_preliminary_H'])
    filtered_e_df = original_e_df[original_e_df['original_hash'].isin(retained_hashes)].copy()

    surviving_p_hashes = set()
    for sources_str in filtered_e_df['sources'].dropna():
        surviving_p_hashes.update(sources_str.split(','))
    surviving_p_hashes = sorted(surviving_p_hashes)
    p_id_map = {p_hash: f"EVODEX.{__version__}-P{idx+1}" for idx, p_hash in enumerate(surviving_p_hashes)}

    filtered_e_df['id'] = [f"EVODEX.{__version__}-E{idx+1}" for idx in range(len(filtered_e_df))]
    # Remove map_e_sources logic; leave sources as raw hashes
    # filtered_e_df['sources'] = filtered_e_df['sources'].apply(map_e_sources)
    filtered_e_df.to_csv(paths['evodex_e_phase3b_final'], index=False)
    print(f"Final pruned EVODEX-E saved to {paths['evodex_e_phase3b_final']}.")

    original_p_df = pd.read_csv(paths['evodex_p_phase3a_retained'])
    filtered_p_df = original_p_df[original_p_df['id'].isin(surviving_p_hashes)].copy()
    # Remove id mapping; keep id as hash
    # filtered_p_df['id'] = filtered_p_df['id'].map(p_id_map)
    filtered_p_df.to_csv(paths['evodex_p_phase3b_final'], index=False)
    print(f"Final pruned EVODEX-P saved to {paths['evodex_p_phase3b_final']}.")

    # Final EVODEX-F
    evodex_f_df = pd.read_csv(paths['evodex_f_filtered'])

    # Define function to check if F row has at least one valid P source
    def f_row_has_valid_source(sources_str):
        if pd.isna(sources_str):
            return False
        sources = sources_str.split(',')
        return any(src in surviving_p_hashes for src in sources)

    # Filter rows where F has at least one valid P source
    filtered_f_df = evodex_f_df[evodex_f_df['sources'].apply(f_row_has_valid_source)].copy()

    # Prune sources to only surviving_p_hashes
    def prune_f_sources(sources_str):
        if pd.isna(sources_str):
            return ''
        sources = sources_str.split(',')
        kept = [src for src in sources if src in surviving_p_hashes]
        return ','.join(kept)

    filtered_f_df['sources'] = filtered_f_df['sources'].apply(prune_f_sources)

    # Save
    filtered_f_df.to_csv(paths['evodex_f_phase3b_final'], index=False)
    print(f"Final pruned EVODEX-F saved to {paths['evodex_f_phase3b_final']}.")

    # Final EVODEX-R
    evodex_r_df = pd.read_csv(paths['evodex_r_preliminary'])
    evodex_r_df['sources'] = evodex_r_df['sources'].astype(str)

    # Gather surviving R hashes from surviving P's sources
    original_p_df = pd.read_csv(paths['evodex_p_phase3a_retained'])
    surviving_p_df = original_p_df[original_p_df['id'].isin(surviving_p_hashes)].copy()
    surviving_r_hashes = set()
    for sources_str in surviving_p_df['sources'].dropna():
        surviving_r_hashes.update(sources_str.split(','))
    surviving_r_hashes = sorted(surviving_r_hashes)

    filtered_r_df = evodex_r_df[evodex_r_df['id'].isin(surviving_r_hashes)].copy()
    filtered_r_df.to_csv(paths['evodex_r_phase3b_final'], index=False)
    print(f"Final pruned EVODEX-R saved to {paths['evodex_r_phase3b_final']}.")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 3b ERO trimming completed in {elapsed_time:.2f} seconds.")


if __name__ == "__main__":
    main()
