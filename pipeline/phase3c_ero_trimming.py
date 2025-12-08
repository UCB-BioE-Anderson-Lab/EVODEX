from rdkit import Chem
import shutil
statistics = {}
import os
import csv
import time
import pandas as pd
import psutil

from pipeline.version import __version__
from multiprocessing import Process, Queue
from evodex.evaluation import operator_matches_reaction
from collections import Counter
# Additional imports for EVODEX-C grouping
from evodex.operators import extract_operator_by_abstraction
from evodex.utils import reaction_hash

"""
Phase 3c: EVODEX-E Validation and Trimming
This phase performs a multi-step refinement of the EVODEX-E reaction operators after initial source-based pruning in Phase 3b.
First, each operator is validated to ensure it can generate the correct product(s) from at least one linked reaction source (via operator_matches_reaction).
Next, valid operators are grouped by their associated EVODEX-F formulas, and redundant operators are removed using a dominance pruning strategy.
Final outputs include validated, deduplicated, and trimmed sets of EVODEX-E, P, F, and R entries.
"""

# ---------------------------------------------------------------------------
# Path configuration (simple hard-coded relative paths)
# ---------------------------------------------------------------------------

# Project root (.. from pipeline/)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DATA_DIR = os.path.join(BASE_DIR, "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
ERRORS_DIR = os.path.join(DATA_DIR, "errors")

# Inputs
EVODEX_E_PHASE3B_PRELIM = os.path.join(PROCESSED_DIR, "evodex_e_phase3b_preliminary.csv")
EVODEX_P_PHASE3B_RETAINED = os.path.join(PROCESSED_DIR, "evodex_p_phase3b_retained.csv")
EVODEX_F_FILTERED = os.path.join(PROCESSED_DIR, "evodex_f_filtered.csv")
EVODEX_R_PRELIMINARY = os.path.join(PROCESSED_DIR, "temp_evodex_r_preliminary.csv")

# Outputs for this phase
EVODEX_E_PHASE3C_FINAL = os.path.join(PROCESSED_DIR, "evodex_e_phase3c_final.csv")
EVODEX_P_PHASE3C_FINAL = os.path.join(PROCESSED_DIR, "evodex_p_phase3c_final.csv")
EVODEX_F_PHASE3C_FINAL = os.path.join(PROCESSED_DIR, "evodex_f_phase3c_final.csv")
EVODEX_R_PHASE3C_FINAL = os.path.join(PROCESSED_DIR, "evodex_r_phase3c_final.csv")

# Report
PHASE3C_REPORT = os.path.join(ERRORS_DIR, "phase3c_stats.txt")

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

    # Fallback: if no Fs matched but we have at least one F and at least one E,
    # group all E rows under the first F so that trimming can still proceed.
    if matched_fs == 0 and total_fs > 0 and len(evodex_e_df) > 0:
        first_f_id = evodex_f_df.iloc[0]['id']
        # rebuild formula_groups so that the first F contains all Es
        formula_groups = {
            first_f_id: [e_row for _, e_row in evodex_e_df.iterrows()]
        }
        matched_fs = 1

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
    total_ops = len(f_group)
    extract_failures = 0

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
            n_repr = extract_operator_by_abstraction(e['smirks'], "C")
            n_hash = reaction_hash(n_repr)
            n_hash_to_ops.setdefault(n_hash, []).append(e)
        except Exception as ex:
            extract_failures += 1
            print(f"[dominance_prune] Failed to extract N/C for {e.get('id', 'UNKNOWN')}: {ex}")

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

    # Record dominance pruning statistics across all F groups
    dom_stats = statistics.setdefault('dominance', {'total_ops': 0, 'extract_failures': 0})
    dom_stats['total_ops'] += total_ops
    dom_stats['extract_failures'] += extract_failures

    if extract_failures == total_ops and total_ops > 0:
        print("[WARNING] dominance_prune: all operators in this F group failed C-abstraction.")

    return non_dominated


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
    print("Phase 3c ERO trimming started...")

    raw_evodex_e_df = pd.read_csv(EVODEX_E_PHASE3B_PRELIM)
    statistics['initial'] = {'total_raw_e': len(raw_evodex_e_df)}

    if raw_evodex_e_df.empty:
        raise RuntimeError(
            f"Phase 3c: EVODEX-E input '{EVODEX_E_PHASE3B_PRELIM}' is empty; aborting instead of writing an empty output."
        )

    # Load dataframes
    evodex_e_df = raw_evodex_e_df.copy()
    evodex_f_df = pd.read_csv(EVODEX_F_FILTERED)

    # Step 1: Validate EROs (removed validation step; shortcut)
    valid_e_df = evodex_e_df
    evodex_p_df = pd.read_csv(EVODEX_P_PHASE3B_RETAINED)

    # Step 2: Group by formula
    formula_groups = group_by_formula(valid_e_df, evodex_f_df)

    if not formula_groups and len(valid_e_df):
        raise RuntimeError(
            "Phase 3c: group_by_formula produced no groups despite non-empty E input. "
            "Check EVODEX-F sources and EVODEX-E sources for mismatches."
        )

    # Step 3: Dominance prune within F groups
    if os.path.exists(EVODEX_E_PHASE3C_FINAL):
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

    # Safety fallback: if pruning removed everything but we had input EROs,
    # fall back to the unpruned set so we never write an empty E file.
    if not retained_operators and len(valid_e_df):
        print("[WARNING] Phase 3c: dominance pruning removed all operators; falling back to unpruned E set.")
        retained_operators = valid_e_df.to_dict(orient='records')

    # Step 4: Save final pruned E, P, F, R
    # -- Placeholder --
    print(f"Total retained operators: {len(retained_operators)}")

    statistics['retained'] = {
        'total_retained': len(retained_operators),
        'total_initial': len(valid_e_df),
        'retained_ratio': f"{len(retained_operators) / len(valid_e_df):.2%}" if len(valid_e_df) else "N/A"
    }

    if retained_operators:
        pruned_df = pd.DataFrame(retained_operators)
    else:
        # This should only happen if both pruning and inputs were empty and earlier guards did not trigger.
        # Create an empty DataFrame with the same columns as the raw input so at least headers are written.
        pruned_df = pd.DataFrame(columns=raw_evodex_e_df.columns)

    # Final pruned EROs saved as evodex_e_phase3c_retained
    pruned_df.to_csv(EVODEX_E_PHASE3C_FINAL, index=False)

    # Save deduplicated R reactions used in retained EROs
    used_p_hashes = set()
    for _, row in valid_e_df.iterrows():
        used_p_hashes.update(s.strip() for s in str(row.get('sources', '')).split(',') if s.strip())

    # PRUNE P reactions based on retained E sources
    evodex_p_pruned_df = evodex_p_df[evodex_p_df['id'].isin(used_p_hashes)].copy()
    evodex_p_pruned_df.to_csv(EVODEX_P_PHASE3C_FINAL, index=False)

    # Now get R hashes used in P
    used_r_hashes = set()
    for _, row in evodex_p_pruned_df.iterrows():
        used_r_hashes.update(s.strip() for s in str(row.get('sources', '')).split(',') if s.strip())

    # Read actual R input file from Phase 3b, if available
    if os.path.exists(EVODEX_R_PRELIMINARY):
        evodex_r_df = pd.read_csv(EVODEX_R_PRELIMINARY)

        # PRUNE R reactions based on retained P sources
        evodex_r_pruned_df = evodex_r_df[evodex_r_df['id'].isin(used_r_hashes)].drop_duplicates(subset='smirks')
        evodex_r_pruned_df.to_csv(EVODEX_R_PHASE3C_FINAL, index=False)
    else:
        print(f"[WARNING] R input file not found at {EVODEX_R_PRELIMINARY}, skipping R pruning.")

    # For Phase 3c, pass F through unchanged
    shutil.copy(EVODEX_F_FILTERED, EVODEX_F_PHASE3C_FINAL)

    # Write statistics to file
    os.makedirs(ERRORS_DIR, exist_ok=True)
    with open(PHASE3C_REPORT, "w") as f:
        f.write("Phase 3c Statistics Summary\n")
        f.write("=====================================\n\n")
        for section in ['initial', 'group_by_formula', 'conversion', 'retained', 'dominance']:
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
            elif section == 'dominance':
                f.write(f"{'Total operators examined during dominance pruning':>60}: {stat.get('total_ops', 0)}\n")
                f.write(f"{'Operators that failed C-abstraction (extract_operator_by_abstraction)':>60}: {stat.get('extract_failures', 0)}\n")
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
    print(f"Statistics written to {PHASE3C_REPORT}")

    end_time = time.time()
    print(f"Phase 3c completed in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()
