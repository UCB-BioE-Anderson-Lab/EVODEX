import shutil
import os
import csv
import pandas as pd
from pipeline.config import load_paths
from evodex.evaluation import match_operators
from evodex.astatine import convert_dataframe_smiles_column

"""
Phase 3b: ERO Trimming

This phase trims the EVODEX-E data to weed out operators from incorrectly mapped reactions
as well as operators with excessive conugation. It does so by using the evaluation.py
algorithms to identify EVODEX-E that satisfy the input EVODEX-P. Based on the table of matches,
it trims out those operators that are dominated by others. Finally, it consolidates the -R, -P, -F, 
and -E datasets to only that data that cleanly results in the non-dominated EVODEX-E.
"""


def main():
    # Load paths
    paths = load_paths('pipeline/config/paths.yaml')
    
    # Load EVODEX-E preliminary into DataFrame
    evodex_e_df = pd.read_csv(paths['evodex_e_preliminary'])
    
    # Convert 'smirks' column At -> H
    converted_df, conversion_errors = convert_dataframe_smiles_column(evodex_e_df, 'smirks')

    # Assign temporary EVODEX-E IDs in the format "EVODEX.1-E<number>_temp"
    temp_ids = [f"EVODEX.1-E{idx+1}_temp" for idx in range(len(converted_df))]
    converted_df['id'] = temp_ids

    # Save converted DataFrame to EVODEX/data
    dst = os.path.join('evodex', 'data', 'EVODEX-E_reaction_operators.csv')
    converted_df.to_csv(dst, index=False)
    print(f"Converted EVODEX-E 'smirks' column and saved to {dst}.")
    
    # Log any conversion errors
    if conversion_errors:
        print(f"Conversion errors encountered in 'smirks' column:")
        for error in conversion_errors:
            print(error)

    # Copy EVODEX-F unique formulas to EVODEX/data
    src_f = paths['evodex_f_filtered']
    dst_f = os.path.join('evodex', 'data', 'EVODEX-F_unique_formulas.csv')
    shutil.copyfile(src_f, dst_f)
    print(f"Copied {src_f} to {dst_f}.")

    # Copy EVODEX-P reaction operators to EVODEX/data
    src_p = paths['evodex_p_filtered']
    dst_p = os.path.join('evodex', 'data', 'EVODEX-P_reaction_operators.csv')
    shutil.copyfile(src_p, dst_p)
    print(f"Copied {src_p} to {dst_p}.")
    
    # Load EVODEX-P filtered data
    evodex_p = pd.read_csv(paths['evodex_p'])

    # Initialize the match table: dict mapping operator IDs to sets of matching EVODEX-P IDs
    match_table = dict()

    # For each EVODEX-P row, run match_operators on its SMIRKS and record matches
    for _, row in evodex_p.iterrows():
        p_id = row['id']
        smirks = row['smirks']
        matched_ops = match_operators(smirks)
        for op_id in matched_ops:
            if op_id not in match_table:
                match_table[op_id] = set()
            match_table[op_id].add(p_id)

    # Save match_table to CSV
    with open('phase3b_match_table.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['operator_id', 'matched_p_ids'])
        for op_id, p_ids in match_table.items():
            writer.writerow([op_id, ';'.join(map(str, sorted(p_ids)))])

    print(f"Phase 3b ERO Trimming scaffold executed.")
    print(f"Total operators matched: {len(match_table)}")
    total_matches = sum(len(pids) for pids in match_table.values())
    print(f"Total matches recorded: {total_matches}")

if __name__ == "__main__":
    main()
