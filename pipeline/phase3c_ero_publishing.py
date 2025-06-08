import csv
import os
import shutil
import time
import pandas as pd
from pipeline.config import load_paths
from pipeline.version import __version__
from evodex.astatine import convert_dataframe_smiles_column
import sys
csv.field_size_limit(sys.maxsize)

"""
Phase 3c: EVODEX Publishing to evodex/data

Goal: Convert final Phase 3b outputs to EVODEX.1 format and publish H-converted versions to evodex/data.

Input:
- evodex_e_phase3b_final (At)
- evodex_p_phase3b_final (At)
- evodex_f_phase3b_final (no SMIRKS, no conversion needed)
- evodex_r_phase3b_final (At)

Steps:
1. Assign final EVODEX.1 IDs to E, P, R.
2. Update source references to match assigned EVODEX-P and EVODEX-R IDs.
3. Write final H-converted E, P, R to:
    - EVODEX-E_reaction_operators.csv
    - EVODEX-P_partial_reactions.csv
    - EVODEX-R_full_reactions.csv
4. Write EVODEX-F_unique_formulas.csv with updated sources.
5. Write raw_data_published (selected_reactions.csv) for website (filtered raw reactions used by EVODEX-R).

Outputs (all in evodex/data/):
- EVODEX-E_reaction_operators.csv
- EVODEX-P_partial_reactions.csv
- EVODEX-R_full_reactions.csv
- EVODEX-F_unique_formulas.csv
- selected_reactions.csv (raw_data_published)
"""

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def assign_ids(df, prefix):
    """Assigns new EVODEX.1 IDs and returns mapping."""
    id_map = {}
    new_rows = []
    for i, row in enumerate(df.itertuples(index=False), start=1):
        new_id = f"EVODEX.{__version__}-{prefix}{i}"
        id_map[row.id] = new_id
        new_row = row._asdict()
        new_row['id'] = new_id
        new_rows.append(new_row)
    return id_map, new_rows

def update_sources(row_sources, source_id_map):
    updated_sources = []
    for src in row_sources.split(','):
        src = src.strip()
        if src in source_id_map:
            updated_sources.append(source_id_map[src])
        else:
            # Optionally log missing source
            pass
    return ','.join(updated_sources)

def main():
    start_time = time.time()
    print("Phase 3c EVODEX publishing started...")
    print("=== Phase 3c: EVODEX Publishing ===")
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # --- Process EVODEX-E ---
    print("Processing EVODEX-E...")
    e_df = pd.read_csv(paths['evodex_e_phase3b_final'])
    e_id_map, e_rows = assign_ids(df=e_df, prefix='E')

    # --- Process EVODEX-P ---
    print("Processing EVODEX-P...")
    p_df = pd.read_csv(paths['evodex_p_phase3b_final'])
    p_id_map, p_rows = assign_ids(df=p_df, prefix='P')

    # Update E sources now that P IDs are known
    e_updated_rows = []
    with open(paths['evodex_e_phase3b_final'], 'r') as efile:
        reader = csv.DictReader(efile)
        for row in reader:
            row['sources'] = update_sources(row['sources'], p_id_map)
            filtered_row = {key: row[key] for key in ['id', 'smirks', 'sources']}
            e_updated_rows.append(filtered_row)

    # --- Process EVODEX-R ---
    print("Processing EVODEX-R...")
    r_df = pd.read_csv(paths['evodex_r_phase3b_final'])
    r_id_map, r_rows = assign_ids(df=r_df, prefix='R')

    # Update P sources now that R IDs are known
    p_updated_rows = []
    with open(paths['evodex_p_phase3b_final'], 'r') as pfile:
        reader = csv.DictReader(pfile)
        for row in reader:
            row['sources'] = update_sources(row['sources'], r_id_map)
            row['id'] = p_id_map[row['id']]
            p_updated_rows.append(row)

    # --- Process EVODEX-F ---
    print("Processing EVODEX-F...")
    f_df = pd.read_csv(paths['evodex_f_phase3b_final'])
    f_id_map, f_rows = assign_ids(df=f_df, prefix='F')

    # Now update F sources (P IDs)
    f_updated_rows = []
    for row in f_rows:
        row['sources'] = update_sources(row['sources'], p_id_map)
        f_updated_rows.append(row)

    # --- Write phase3c_final files and publish ---

    # EVODEX-E: Write updated rows as-is to phase3c_final and publish unchanged
    print("Writing EVODEX-E phase3c_final file...")
    with open(paths['evodex_e_phase3c_final'], 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(e_updated_rows)
    print("Publishing EVODEX-E to evodex/data...")
    shutil.copyfile(paths['evodex_e_phase3c_final'], os.path.join('evodex/data/EVODEX-E_reaction_operators.csv'))

    # EVODEX-P: Write updated rows to phase3c_final, convert to H for publishing
    print("Writing EVODEX-P phase3c_final file...")
    with open(paths['evodex_p_phase3c_final'], 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(p_updated_rows)
    print("Publishing EVODEX-P to evodex/data...")
    p_df_h, _ = convert_dataframe_smiles_column(pd.DataFrame(p_updated_rows), 'smirks')
    with open(os.path.join('evodex/data/EVODEX-P_partial_reactions.csv'), 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'smirks', 'sources'])
        writer.writeheader()
        writer.writerows(p_df_h.to_dict(orient='records'))

    # EVODEX-R: Write updated rows to phase3c_final, convert to H for publishing
    print("Writing EVODEX-R phase3c_final file...")
    with open(paths['evodex_r_phase3c_final'], 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(r_rows)
    print("Publishing EVODEX-R to evodex/data...")
    r_df_h, _ = convert_dataframe_smiles_column(pd.DataFrame(r_rows), 'smirks')
    with open(os.path.join('evodex/data/EVODEX-R_full_reactions.csv'), 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'smirks', 'sources'])
        writer.writeheader()
        writer.writerows(r_df_h.to_dict(orient='records'))

    # EVODEX-F: Write updated rows to phase3c_final and publish unchanged
    print("Writing EVODEX-F phase3c_final file...")
    with open(paths['evodex_f_phase3c_final'], 'w', newline='') as outfile:
        reader_fieldnames = f_updated_rows[0].keys() if f_updated_rows else ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=reader_fieldnames)
        writer.writeheader()
        writer.writerows(f_updated_rows)
    print("Publishing EVODEX-F to evodex/data...")
    shutil.copyfile(paths['evodex_f_phase3c_final'], os.path.join('evodex/data/EVODEX-F_unique_formulas.csv'))

    # --- Write raw_data_published ---
    print("Writing raw_data_published (selected_reactions.csv)...")

    # Load EVODEX-R sources
    r_df_final = pd.read_csv(paths['evodex_r_phase3c_final'], dtype={'sources': str})
    source_ids = set()
    for sources_str in r_df_final['sources']:
        for src in sources_str.split(','):
            src = src.strip()
            if src:
                source_ids.add(src)

    # Load raw_data (not filtered_data), read rxn_idx as string
    raw_data_df = pd.read_csv(paths['raw_data'], dtype={'rxn_idx': str})

    # Filter rows where rxn_idx is in source_ids
    raw_data_published_df = raw_data_df[raw_data_df['rxn_idx'].isin(source_ids)].copy()
    # De-duplicate by rxn_idx so only one row per rxn_idx is kept
    raw_data_published_df = raw_data_published_df.drop_duplicates(subset='rxn_idx')

    # Write to raw_data_published path (keep same columns)
    raw_data_published_df.to_csv(paths['raw_data_published'], index=False)

    print(f"Published raw_data_published to {paths['raw_data_published']}")

    print("=== Phase 3c publishing complete ===")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 3c EVODEX publishing completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
