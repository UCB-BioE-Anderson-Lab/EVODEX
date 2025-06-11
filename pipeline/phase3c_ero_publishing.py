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

Goal: Convert final Phase 3b outputs to EVODEX.1 format and publish cleaned, H-converted versions to evodex/data.

Assumes the following dependency structure:
  - EVODEX-R is the foundation (source for P)
  - EVODEX-P uses R as sources
  - EVODEX-E and F use P as sources

Processing steps:
  1. Assign final EVODEX.1 IDs to R → P → E → F in dependency order.
  2. Replace source hashes using previously assigned ID maps (R→P, P→E/F).
  3. Publish each EVODEX file with updated IDs and sources.
  4. Apply At→H conversion for EVODEX-P and R SMIRKS fields.
  5. Publish raw_data subset used by EVODEX-R.
  6. Cleanup legacy EVODEX data JSONs.
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

def update_sources(row_sources, source_id_map, row_id=None):
    updated_sources = []
    for src in row_sources.split(','):
        src = src.strip()
        if src in source_id_map:
            updated_sources.append(source_id_map[src])
        else:
            print(f"[ERROR] Missing source mapping in row {row_id or 'UNKNOWN'}: '{src}' not found in source_id_map")
            raise ValueError(f"Unresolved source '{src}' in row {row_id or 'UNKNOWN'}")
    return ','.join(updated_sources)

def main():
    start_time = time.time()
    print("Phase 3c EVODEX publishing started...")
    print("=== Phase 3c: EVODEX Publishing ===")
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # --- Process EVODEX-R (first) ---
    print("Processing EVODEX-R...")
    r_df = pd.read_csv(paths['evodex_r_phase3b_final'])
    r_id_map, r_rows = assign_ids(df=r_df, prefix='R')

    # --- Process EVODEX-P (second) ---
    print("Processing EVODEX-P...")
    p_df = pd.read_csv(paths['evodex_p_phase3b_final'])
    p_id_map, p_rows = assign_ids(df=p_df, prefix='P')

    # Update P sources now that R IDs are known
    p_updated_rows = []
    with open(paths['evodex_p_phase3b_final'], 'r') as pfile:
        reader = csv.DictReader(pfile)
        for row in reader:
            row['sources'] = update_sources(row['sources'], r_id_map, row['id'])
            row['id'] = p_id_map[row['id']]
            p_updated_rows.append(row)

    # --- Process EVODEX-E (third) ---
    print("Processing EVODEX-E...")
    e_df = pd.read_csv(paths['evodex_e_phase3b_final'])
    e_id_map, e_rows = assign_ids(df=e_df, prefix='E')

    # Update E sources now that P IDs are known
    e_updated_rows = []
    with open(paths['evodex_e_phase3b_final'], 'r') as efile:
        reader = csv.DictReader(efile)
        for row in reader:
            original_id = row['id']
            row['id'] = e_id_map[original_id]
            row['sources'] = update_sources(row['sources'], p_id_map, row['id'])
            filtered_row = {key: row[key] for key in ['id', 'smirks', 'sources']}
            e_updated_rows.append(filtered_row)

    # --- Process EVODEX-F (fourth) ---
    print("Processing EVODEX-F...")
    f_df = pd.read_csv(paths['evodex_f_phase3b_final'])
    f_updated_rows = []
    for i, row in enumerate(f_df.itertuples(index=False), start=1):
        new_id = f"EVODEX.{__version__}-F{i}"
        updated_sources = update_sources(row.sources, p_id_map, new_id)
        f_updated_rows.append({
            'id': new_id,
            'formula': row.formula,
            'sources': updated_sources
        })

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
        reader_fieldnames = f_updated_rows[0].keys() if f_updated_rows else ['id', 'formula', 'sources']
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

    # --- Cleanup old EVODEX data JSON files ---
    json_files = [
        "evodex/data/evaluation_operator_data.json",
        "evodex/data/evodex_e_data.json"
    ]
    for file in json_files:
        if os.path.exists(file):
            os.remove(file)
    print("Cleaned up old EVODEX data JSON files.")

    print("=== Phase 3c publishing complete ===")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 3c EVODEX publishing completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
