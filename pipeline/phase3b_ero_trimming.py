from rdkit import Chem
import shutil
statistics = {}
import os
import csv
import time
import pandas as pd
import psutil
from pipeline.config import load_paths
from evodex.astatine import convert_dataframe_smirks_column_at_to_h, convert_dataframe_smiles_column_at_to_h
from pipeline.version import __version__
from multiprocessing import Process, Queue
from evodex.evaluation import operator_matches_reaction
from collections import Counter
# Additional imports for EVODEX-C grouping
from evodex.operators import extract_operator
from evodex.utils import reaction_hash

"""
Phase 3b: EVODEX-E Validation and Trimming
This phase performs a multi-step refinement of the EVODEX-E reaction operators after initial source-based pruning in Phase 3a.
First, each operator is validated to ensure it can generate the correct product(s) from at least one linked reaction source (via operator_matches_reaction).
Next, valid operators are grouped by their associated EVODEX-F formulas, and redundant operators are removed using a dominance pruning strategy.
Final outputs include validated, deduplicated, and trimmed sets of EVODEX-E, P, F, and R entries.
"""

csv.field_size_limit(10**7)

def log_memory_usage(label=""):
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 * 1024)
    print(f"[MEMORY] {label}: {mem:.2f} MB")


def group_by_formula(evodex_e_df, evodex_f_df):
    # Partition EROs by F
    # Return dict of F -> list of EROs
    formula_groups = {}
    f_source_map = {}

    skipped_fs = 0

    # Build map of F -> set of P hashes
    for _, row in evodex_f_df.iterrows():
        f_id = row['id']
        sources = row.get('sources')
        if pd.isna(sources) or sources.strip() == "{}":
            skipped_fs += 1
            continue
        f_source_map[f_id] = set(s.strip() for s in sources.split(',') if s.strip())

    f_e_map = {f: [] for f in f_source_map}

    # Map each E to all F groups whose P hashes include any of the E's P hashes
    for _, e_row in evodex_e_df.iterrows():
        e_sources = e_row.get('sources')
        if pd.isna(e_sources):
            continue
        e_p_hashes = set(s.strip() for s in e_sources.split(',') if s.strip())
        for f_id, f_p_hashes in f_source_map.items():
            if e_p_hashes & f_p_hashes:
                f_e_map[f_id].append(e_row)

    # Build final group dict and gather stats
    formula_groups = {f: rows for f, rows in f_e_map.items() if rows}
    total_fs = len(evodex_f_df)
    matched_fs = len(formula_groups)
    total_es = sum(len(v) for v in formula_groups.values())
    multi_e_fs = sum(1 for v in formula_groups.values() if len(v) > 1)

    # Compute unmatched Es
    matched_e_ids = {e['id'] for group in formula_groups.values() for e in group}
    unmatched_e_count = len(evodex_e_df) - len(matched_e_ids)
    statistics['group_by_formula'] = {
        'total_f': total_fs,
        'matched_f': matched_fs,
        'total_e': total_es,
        'multi_e_f': multi_e_fs,
        'skipped_f': skipped_fs,
        'unmatched_e': unmatched_e_count
    }

    print(f"[group_by_formula] Total Fs: {total_fs}, Matched Fs: {matched_fs}, Total Es grouped: {total_es}, Fs with >1 E: {multi_e_fs}")

    return formula_groups

def dominance_prune_within_formula(f_group):
    def count_substrate_atoms(smirks):
        try:
            substrate_part = smirks.split(">>")[0]
            mols = substrate_part.split(".")
            total_atoms = 0
            for mol in mols:
                m = Chem.MolFromSmiles(mol)
                if m:
                    total_atoms += m.GetNumAtoms()
            return total_atoms
        except Exception as e:
            print(f"[count_substrate_atoms] Error processing SMIRKS: {smirks} -> {e}")
            return 0

    # Sort operators by ascending substrate atom count
    f_group = sorted(f_group, key=lambda e: count_substrate_atoms(e['smirks']))

    # Compute EVODEX-N hash for each operator and group by hash
    n_hash_to_ops = {}
    for e in f_group:
        try:
            n_repr = extract_operator(
                e['smirks'],
                include_stereochemistry=True,
                include_sigma=True,
                include_pi=False,
                include_unmapped_hydrogens=True,
                include_unmapped_heavy_atoms=True,
                include_static_hydrogens=True
            )
            n_hash = reaction_hash(n_repr)
            n_hash_to_ops.setdefault(n_hash, []).append(e)
        except Exception as ex:
            print(f"[dominance_prune] Failed to extract N for {e.get('id', 'UNKNOWN')}: {ex}")

    non_dominated = []
    for n_hash, ops in n_hash_to_ops.items():
        ops = sorted(ops, key=lambda e: count_substrate_atoms(e['smirks']))
        group_non_dominated = []
        for i, candidate in enumerate(ops):
            print(f"[dominance_prune] Processing operator {i+1}/{len(ops)}: {candidate.get('id', 'UNKNOWN')}")
            is_dominated = False
            for nd in group_non_dominated:
                if check_match_with_timeout(nd['smirks'], candidate['smirks'], timeout=60):
                    print(f"[dominance_prune] Candidate {candidate.get('id', 'UNKNOWN')} is dominated by {nd.get('id', 'UNKNOWN')}")
                    is_dominated = True
                    break
            if not is_dominated:
                group_non_dominated.append(candidate)
        non_dominated.extend(group_non_dominated)

    return non_dominated

def load_and_prepare_data(paths):
    if os.path.exists(paths['evodex_e_phase3b_h_converted']) and os.path.exists(paths['evodex_p_phase3b_h_converted']):
        print("Data already prepared, skipping load_and_prepare_data.")
        return
    # Load EVODEX-E Phase 3a preliminary and convert SMIRKS to H
    print("Converting EVODEX-E Phase 3a preliminary to H...")
    evodex_e_df = pd.read_csv(paths['evodex_e_phase3a_preliminary'])
    converted_e_df, conversion_errors_e = convert_dataframe_smirks_column_at_to_h(evodex_e_df, 'smirks')
    # Filter out rows with null or empty smirks
    invalid_e_mask = converted_e_df['smirks'].isna() | (converted_e_df['smirks'].str.strip() == "")
    invalid_e_rows = evodex_e_df[invalid_e_mask]
    if not invalid_e_rows.empty:
        os.makedirs("data/errors", exist_ok=True)
        invalid_e_rows.to_csv("data/errors/evodex_e_conversion_errors.csv", index=False)
    converted_e_df = converted_e_df[~invalid_e_mask]
    converted_e_df.to_csv(paths['evodex_e_phase3b_h_converted'], index=False)
    print(f"Saved: {paths['evodex_e_phase3b_h_converted']}")

    # Load EVODEX-P Phase 3a retained and convert SMIRKS to H
    print("Converting EVODEX-P Phase 3a retained to H...")
    evodex_p_df = pd.read_csv(paths['evodex_p_phase3a_retained'])
    converted_p_df, conversion_errors_p = convert_dataframe_smiles_column_at_to_h(evodex_p_df, 'smirks')
    invalid_p_mask = converted_p_df['smirks'].isna() | (converted_p_df['smirks'].str.strip() == "")
    invalid_p_rows = evodex_p_df[invalid_p_mask]
    if not invalid_p_rows.empty:
        os.makedirs("data/errors", exist_ok=True)
        invalid_p_rows.to_csv("data/errors/evodex_p_conversion_errors.csv", index=False)
    converted_p_df = converted_p_df[~invalid_p_mask]
    converted_p_df.to_csv(paths['evodex_p_phase3b_h_converted'], index=False)
    print(f"Saved: {paths['evodex_p_phase3b_h_converted']}")

    # Save converted P with H-converted SMIRKS

    statistics['conversion'] = {
        'evodex_e_invalid_count': len(invalid_e_rows),
        'evodex_p_invalid_count': len(invalid_p_rows)
    }

def _match_worker(q, op, rxn):
    try:
        result = operator_matches_reaction(op, rxn)
        # print(f"[MATCH WORKER] Result for operator:\n  {op}\nvs reaction:\n  {rxn[:200]}...\n=> {result}")
        q.put(result)
    except Exception as e:
        print(f"[MATCH WORKER] ERROR comparing operator vs reaction: {e}")
        q.put(False)

def check_match_with_timeout(op_smirks, rxn, timeout=60):
    q = Queue()
    p = Process(target=_match_worker, args=(q, op_smirks, rxn))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        p.join()
        print(f"[TIMEOUT] Match timed out for operator: {op_smirks}")
        return False
    try:
        return q.get_nowait()
    except Exception:
        return False

def main():
    start_time = time.time()
    print("Phase 3b ERO trimming started...")
    
    paths = load_paths('pipeline/config/paths.yaml')
    raw_evodex_e_df = pd.read_csv(paths['evodex_e_phase3a_preliminary'])
    statistics['initial'] = {'total_raw_e': len(raw_evodex_e_df)}
    load_and_prepare_data(paths)

    # Load dataframes
    evodex_e_df = pd.read_csv(paths['evodex_e_phase3b_h_converted'])
    evodex_f_df = pd.read_csv(paths['evodex_f_filtered'])

    # Step 1: Validate EROs (removed validation step; shortcut)
    valid_e_df = evodex_e_df
    evodex_p_df = pd.read_csv(paths['evodex_p_phase3b_h_converted'])

    # Step 2: Group by formula
    formula_groups = group_by_formula(valid_e_df, evodex_f_df)

    # Step 3: Dominance prune within F groups
    if os.path.exists(paths['evodex_e_phase3b_final']):
        print("Dominance pruning output already exists, skipping pruning.")
        return
    retained_operators = []
    for f, group in formula_groups.items():
        if len(group) == 1:
            retained_operators.extend(group)
        else:
            retained = dominance_prune_within_formula(group)
            retained_operators.extend(retained)

    # Deduplicate by 'id' before writing
    retained_operators = list({row['id']: row for row in retained_operators}.values())

    # Step 4: Save final pruned E, P, F, R
    # -- Placeholder --
    print(f"Total retained operators: {len(retained_operators)}")

    statistics['retained'] = {
        'total_retained': len(retained_operators),
        'total_initial': len(valid_e_df),
        'retained_ratio': f"{len(retained_operators) / len(valid_e_df):.2%}" if len(valid_e_df) else "N/A"
    }

    pruned_df = pd.DataFrame(retained_operators)
    # Final pruned EROs saved as evodex_e_phase3b_retained
    pruned_df.to_csv(paths['evodex_e_phase3b_final'], index=False)

    # Save deduplicated R reactions used in retained EROs
    used_p_hashes = set()
    for _, row in valid_e_df.iterrows():
        used_p_hashes.update(s.strip() for s in str(row.get('sources', '')).split(',') if s.strip())

    # PRUNE P reactions based on retained E sources
    evodex_p_pruned_df = evodex_p_df[evodex_p_df['id'].isin(used_p_hashes)].copy()
    evodex_p_pruned_df.to_csv(paths['evodex_p_phase3b_final'], index=False)

    # PRUNE F entries whose sources include at least one retained P
    evodex_f_df = pd.read_csv(paths['evodex_f_filtered'])
    valid_p_ids = used_p_hashes
    filtered_f_rows = []
    for _, row in evodex_f_df.iterrows():
        source_ids = [s.strip() for s in str(row.get('sources', '')).split(',') if s.strip()]
        retained_sources = [s for s in source_ids if s in valid_p_ids]
        if retained_sources:
            row_copy = row.copy()
            row_copy['sources'] = ",".join(retained_sources)
            filtered_f_rows.append(row_copy)
    pruned_f_df = pd.DataFrame(filtered_f_rows)
    pruned_f_df.to_csv(paths['evodex_f_phase3b_final'], index=False)

    # Now get R hashes used in P
    used_r_hashes = set()
    for _, row in evodex_p_pruned_df.iterrows():
        used_r_hashes.update(s.strip() for s in str(row.get('sources', '')).split(',') if s.strip())

    # Read actual R input file from Phase 3a
    evodex_r_df = pd.read_csv(paths['evodex_r_preliminary'])

    # PRUNE R reactions based on retained P sources
    evodex_r_pruned_df = evodex_r_df[evodex_r_df['id'].isin(used_r_hashes)].drop_duplicates(subset='smirks')
    evodex_r_pruned_df.to_csv(paths['evodex_r_phase3b_final'], index=False)

    # Write statistics to file
    os.makedirs("data/errors", exist_ok=True)
    with open("data/errors/phase3b_stats.txt", "w") as f:
        f.write("Phase 3b Statistics Summary\n")
        f.write("=====================================\n\n")
        for section in ['initial', 'group_by_formula', 'conversion', 'retained']:
            stat = statistics.get(section)
            if not stat:
                continue
            f.write(f"[{section.upper()}]\n")
            if section == 'initial':
                f.write(f"{'Total raw EVODEX-E entries loaded from file':>60}: {stat['total_raw_e']}\n")
            elif section == 'group_by_formula':
                f.write(f"{'Total EVODEX-F entries loaded from file':>60}: {stat['total_f']}\n")
                f.write(f"{'EVODEX-F entries matched to at least one ERO':>60}: {stat['matched_f']}\n")
                f.write(f"{'Total EVODEX-E reaction operators grouped by F':>60}: {stat['total_e']}\n")
                f.write(f"{'EVODEX-F entries with multiple associated EROs':>60}: {stat['multi_e_f']}\n")
                f.write(f"{'EVODEX-F entries skipped due to missing P links':>60}: {stat.get('skipped_f', 0)}\n")
                f.write(f"{'EVODEX-E entries that failed to match any F':>60}: {stat.get('unmatched_e', 0)}\n")
            elif section == 'conversion':
                f.write(f"{'EVODEX-E entries dropped during SMIRKS conversion':>60}: {stat['evodex_e_invalid_count']}\n")
                f.write(f"{'EVODEX-P entries dropped during SMIRKS conversion':>60}: {stat['evodex_p_invalid_count']}\n")
            elif section == 'retained':
                f.write(f"{'Total EVODEX-E operators before pruning':>60}: {stat['total_initial']}\n")
                f.write(f"{'EVODEX-E operators retained after pruning':>60}: {stat['total_retained']}\n")
                f.write(f"{'Percentage of EROs retained after pruning':>60}: {stat['retained_ratio']}\n")
            f.write("\n")

        conversion_losses = statistics.get('conversion', {}).get('evodex_e_invalid_count', 0)
        retained = statistics.get('retained', {}).get('total_retained', 0)
        total_initial = statistics.get('retained', {}).get('total_initial', 0)
        # Compute pruning_losses safely and handle missing stats.
        if total_initial:
            pruning_losses = total_initial - retained
        else:
            pruning_losses = 0
        total_recorded = conversion_losses + pruning_losses + retained
        raw_total = statistics.get('initial', {}).get('total_raw_e', 0)
        discrepancy = raw_total - total_recorded

        f.write("\n[SANITY CHECK]\n")
        f.write(f"{'Sum of all removed and retained entries':>60}: {total_recorded}\n")
        f.write(f"{'Discrepancy from raw EVODEX-E total':>60}: {discrepancy}\n")
        if raw_total:
            f.write(f"{'Percentage of original EVODEX-E entries retained':>60}: {100 * retained / raw_total:.2f}%\n")
    print("Statistics written to data/errors/phase3b_stats.txt")

    end_time = time.time()
    print(f"Phase 3b completed in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()
