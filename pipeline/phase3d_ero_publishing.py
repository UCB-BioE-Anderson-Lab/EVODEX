import csv
import os
import shutil
import time
import pandas as pd
from pipeline.version import __version__
from evodex.astatine import  convert_dataframe_smiles_column_at_to_h
import sys
csv.field_size_limit(sys.maxsize)

# ---------------------------------------------------------------------------
# Path configuration (simple hard-coded relative paths)
# ---------------------------------------------------------------------------

# Project root (.. from pipeline/)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
RAW_DATA_DIR = os.path.join(DATA_DIR, "raw")
EVODEX_DATA_DIR = os.path.join(BASE_DIR, "evodex", "data")

# Phase 3b final inputs
EVODEX_E_PHASE3C_FINAL = os.path.join(PROCESSED_DIR, "evodex_e_phase3c_final.csv")
EVODEX_P_PHASE3C_FINAL = os.path.join(PROCESSED_DIR, "evodex_p_phase3c_final.csv")
EVODEX_F_PHASE3C_FINAL = os.path.join(PROCESSED_DIR, "evodex_f_phase3c_final.csv")
EVODEX_R_PHASE3C_FINAL = os.path.join(PROCESSED_DIR, "evodex_r_phase3c_final.csv")

# Phase 3d processed outputs (At versions, for future pipeline steps)
EVODEX_P_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_p_phase3d_final.csv")
EVODEX_R_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_r_phase3d_final.csv")
EVODEX_F_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_f_phase3d_final.csv")
EVODEX_E_PHASE3D_FINAL = os.path.join(PROCESSED_DIR, "evodex_e_phase3d_final.csv")

# Raw data and published subset
RAW_DATA = os.path.join(RAW_DATA_DIR, "raw_reactions.csv")
RAW_DATA_PUBLISHED = os.path.join(EVODEX_DATA_DIR, "selected_reactions.csv")

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

def assign_ids(df, prefix):
    """Assigns new EVODEX.1 IDs and returns mapping."""
    id_map = {}
    new_rows = []
    for i, row in enumerate(df.itertuples(index=False), start=1):
        new_id = f"EVODEX.{__version__}-{prefix}{i}"
        id_map[row.id] = new_id
        new_row = row._asdict()
        new_row['id'] = new_id
        new_row['original_id'] = row.id
        new_rows.append(new_row)
    return id_map, new_rows

def update_sources(row_sources, source_id_map, row_id=None, allow_missing=False):
    """
    Map source hashes using source_id_map.

    If allow_missing is False (default), any missing source raises an error.
    If allow_missing is True, missing sources are silently dropped; if no
    sources remain after mapping, None is returned.
    """
    updated_sources = []
    for src in row_sources.split(','):
        src = src.strip()
        if not src:
            continue
        if src in source_id_map:
            updated_sources.append(source_id_map[src])
        else:
            if allow_missing:
                # For F: we may have sources from pruned P; just skip them.
                continue
            print(f"[ERROR] Missing source mapping in row {row_id or 'UNKNOWN'}: '{src}' not found in source_id_map")
            raise ValueError(f"Unresolved source '{src}' in row {row_id or 'UNKNOWN'}")
    if not updated_sources and allow_missing:
        # Caller can decide to drop this row entirely.
        print(f"[WARNING] All sources for row {row_id or 'UNKNOWN'} were missing; dropping this operator.")
        return None
    return ','.join(updated_sources)

def main():
    start_time = time.time()
    print("Phase 3d EVODEX publishing started...")
    print("=== Phase 3d: EVODEX Publishing ===")

    # --- Process EVODEX-R (first) ---
    print("Processing EVODEX-R...")
    r_df = pd.read_csv(EVODEX_R_PHASE3C_FINAL)
    r_id_map, r_rows = assign_ids(df=r_df, prefix='R')

    # --- Process EVODEX-P (second) ---
    print("Processing EVODEX-P...")
    p_df = pd.read_csv(EVODEX_P_PHASE3C_FINAL)
    p_id_map, p_rows = assign_ids(df=p_df, prefix='P')

    # Update P sources now that R IDs are known
    p_updated_rows = []
    with open(EVODEX_P_PHASE3C_FINAL, 'r') as pfile:
        reader = csv.DictReader(pfile)
        for row in reader:
            original_id = row['id']
            row['original_id'] = original_id
            row['sources'] = update_sources(row['sources'], r_id_map, original_id)
            row['id'] = p_id_map[original_id]
            p_updated_rows.append(row)

    # --- Process EVODEX-E (third) ---
    print("Processing EVODEX-E...")
    e_df = pd.read_csv(EVODEX_E_PHASE3C_FINAL)
    e_id_map, e_rows = assign_ids(df=e_df, prefix='E')

    # Update E sources now that P IDs are known
    e_updated_rows = []
    with open(EVODEX_E_PHASE3C_FINAL, 'r') as efile:
        reader = csv.DictReader(efile)
        for row in reader:
            original_id = row['id']
            row['id'] = e_id_map[original_id]
            row['original_id'] = original_id
            row['sources'] = update_sources(row['sources'], p_id_map, row['id'])
            filtered_row = {key: row[key] for key in ['id', 'original_id', 'smirks', 'sources']}
            e_updated_rows.append(filtered_row)

    # --- Process EVODEX-F (fourth) ---
    print("Processing EVODEX-F...")
    f_df = pd.read_csv(EVODEX_F_PHASE3C_FINAL)
    f_updated_rows = []
    for i, row in enumerate(f_df.itertuples(index=False), start=1):
        new_id = f"EVODEX.{__version__}-F{i}"
        updated_sources = update_sources(row.sources, p_id_map, new_id, allow_missing=True)
        if updated_sources is None:
            # All sources came from pruned P; drop this F operator.
            continue
        f_updated_rows.append({
            'id': new_id,
            'original_id': row.id,
            'formula_diff': row.formula_diff,
            'sources': updated_sources
        })

    # --- Write phase3d_final files and publish ---

    # EVODEX-E: Write updated rows as-is to phase3d_final and publish unchanged
    print("Writing EVODEX-E phase3d_final file...")
    with open(EVODEX_E_PHASE3D_FINAL, 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in e_updated_rows:
            writer.writerow({key: row[key] for key in fieldnames})
    print("Publishing EVODEX-E to evodex/data...")
    shutil.copyfile(EVODEX_E_PHASE3D_FINAL, os.path.join(EVODEX_DATA_DIR, "EVODEX-E_reaction_operators.csv"))

    # EVODEX-P: Write updated rows to phase3d_final, convert to H for publishing
    print("Writing EVODEX-P phase3d_final file...")
    with open(EVODEX_P_PHASE3D_FINAL, 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'formula_diff', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in p_updated_rows:
            writer.writerow({key: row.get(key, '') for key in fieldnames})
    print("Publishing EVODEX-P to evodex/data...")
    p_df_h, _ = convert_dataframe_smiles_column_at_to_h(pd.DataFrame(p_updated_rows), 'smirks')
    # Drop original_id from published P
    for col in list(p_df_h.columns):
        if col not in ['id', 'smirks', 'formula_diff', 'sources']:
            p_df_h = p_df_h.drop(columns=[col])
    with open(os.path.join(EVODEX_DATA_DIR, "EVODEX-P_partial_reactions.csv"), 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'formula_diff', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(p_df_h.to_dict(orient='records'))

    # EVODEX-R: Write updated rows to phase3d_final, convert to H for publishing
    print("Writing EVODEX-R phase3d_final file...")
    with open(EVODEX_R_PHASE3D_FINAL, 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in r_rows:
            writer.writerow({key: row[key] for key in fieldnames})
    print("Publishing EVODEX-R to evodex/data...")
    r_df_h, _ = convert_dataframe_smiles_column_at_to_h(pd.DataFrame(r_rows), 'smirks')
    # Drop original_id from published R
    for col in list(r_df_h.columns):
        if col not in ['id', 'smirks', 'sources']:
            r_df_h = r_df_h.drop(columns=[col])
    with open(os.path.join(EVODEX_DATA_DIR, "EVODEX-R_full_reactions.csv"), 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(r_df_h.to_dict(orient='records'))

    # EVODEX-F: Write updated rows to phase3d_final and publish unchanged
    print("Writing EVODEX-F phase3d_final file...")
    with open(EVODEX_F_PHASE3D_FINAL, 'w', newline='') as outfile:
        fieldnames = ['id', 'formula_diff', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in f_updated_rows:
            writer.writerow({key: row[key] for key in fieldnames})
    print("Publishing EVODEX-F to evodex/data...")
    shutil.copyfile(EVODEX_F_PHASE3D_FINAL, os.path.join(EVODEX_DATA_DIR, "EVODEX-F_unique_formulas.csv"))

    # --- Write raw_data_published ---
    print("Writing raw_data_published (selected_reactions.csv)...")

    # Load EVODEX-R sources
    r_df_final = pd.read_csv(EVODEX_R_PHASE3D_FINAL, dtype={'sources': str})
    source_ids = set()
    for sources_str in r_df_final['sources']:
        for src in sources_str.split(','):
            src = src.strip()
            if src:
                source_ids.add(src)

    # Load raw_data (not filtered_data), read rxn_idx as string
    raw_data_df = pd.read_csv(RAW_DATA, dtype={'rxn_idx': str})

    # Filter rows where rxn_idx is in source_ids
    raw_data_published_df = raw_data_df[raw_data_df['rxn_idx'].isin(source_ids)].copy()
    # De-duplicate by rxn_idx so only one row per rxn_idx is kept
    raw_data_published_df = raw_data_published_df.drop_duplicates(subset='rxn_idx')

    # Write to raw_data_published path (keep same columns)
    raw_data_published_df.to_csv(RAW_DATA_PUBLISHED, index=False)

    print(f"Published raw_data_published to {RAW_DATA_PUBLISHED}")

    # --- Cleanup old EVODEX data JSON files ---
    json_files = [
        os.path.join(EVODEX_DATA_DIR, "evaluation_operator_data.json"),
        os.path.join(EVODEX_DATA_DIR, "evodex_e_data.json"),
    ]
    for file in json_files:
        if os.path.exists(file):
            os.remove(file)
    print("Cleaned up old EVODEX data JSON files.")

    print("=== Phase 3d publishing complete ===")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 3d EVODEX publishing completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
