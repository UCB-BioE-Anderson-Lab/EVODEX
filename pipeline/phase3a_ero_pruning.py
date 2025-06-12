import csv
import os
import time
from collections import defaultdict
from pipeline.config import load_paths
from pipeline.version import __version__
import sys
csv.field_size_limit(sys.maxsize)
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Phase 3a: EVODEX-E Pruning
# This phase prunes EVODEX-E operators from the full EVODEX-E set.
# Operators whose reactants lack any atom-mapped atoms are excluded.
# Among the rest, operators with >=5 sources are retained; others are excluded.
# Up to 5 P's are kept per retained E.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    start_time = time.time()
    print("Phase 3a EVODEX-E pruning started...")
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # Load EVODEX-E full entries
    e_entries = []
    with open(paths['evodex_e_phase3_all'], 'r') as efile:
        reader = csv.DictReader(efile)
        for row in reader:
            e_entries.append(row)

    # Build operator_map from EVODEX-E entries
    operator_map = defaultdict(lambda: {'smirks': None, 'sources': set()})
    for row in e_entries:
        op_hash = row['id']
        operator_map[op_hash]['smirks'] = row['smirks']
        if row['sources']:
            sources_set = set(row['sources'].split(','))
            operator_map[op_hash]['sources'].update(sources_set)

    def has_mapped_atoms(op_id, smirks):
        try:
            rxn = AllChem.ReactionFromSmarts(smirks)
            for mol in rxn.GetReactants():
                if mol is None:
                    continue
                for atom in mol.GetAtoms():
                    if atom.GetAtomMapNum() > 0:
                        return True
            print(f"Rejected operator {op_id} — no mapped atoms found.")
            return False
        except Exception as e:
            print(f"Error parsing SMIRKS for operator {op_id}: {e}")
            return False

    rejected_unmapped = 0
    filtered_operator_map = {}
    for k, v in operator_map.items():
        if v['smirks'] and has_mapped_atoms(k, v['smirks']):
            filtered_operator_map[k] = v
        else:
            rejected_unmapped += 1
    operator_map = filtered_operator_map

    # Retain operators with >=5 sources (operators with fewer than 5 sources are excluded)
    retained_ops = {k: v for k, v in operator_map.items() if len(v['sources']) >= 5}

    # Statistics on retained operators
    import statistics
    num_sources_list = [len(data['sources']) for data in retained_ops.values()]
    min_sources = min(num_sources_list) if num_sources_list else 0
    max_sources = max(num_sources_list) if num_sources_list else 0
    mean_sources = statistics.mean(num_sources_list) if num_sources_list else 0
    median_sources = statistics.median(num_sources_list) if num_sources_list else 0

    # Retention breakdown (dropping operators with only 1 source, keeping up to 5 sources per retained E)
    num_ops_excluded = sum(1 for v in operator_map.values() if len(v['sources']) <= 5)
    num_ops_retained = len(retained_ops)
    num_ops_trimmed = sum(1 for v in retained_ops.values() if len(v['sources']) > 5)

    print("Operator retention breakdown:")
    print(f"  Operators with <5 sources (excluded): {num_ops_excluded}")
    print(f"  Operators retained (>=5 sources): {num_ops_retained}")
    print(f"  Operators with >5 sources and thus trimmed to 5: {num_ops_trimmed}")
    print(f"  Operators rejected for missing mapped atoms: {rejected_unmapped}")

    # Sort retained_ops by number of sources (descending)
    sorted_retained_ops = sorted(
        retained_ops.items(),
        key=lambda item: len(item[1]['sources']),
        reverse=True
    )

    # Write retained EVODEX-E operators to preliminary file (trim sources to up to 5 per operator)
    with open(paths['evodex_e_phase3a_preliminary'], 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for op_hash, data in sorted_retained_ops:
            selected_sources = sorted(data['sources'])[:5] if data['sources'] else []
            writer.writerow({
                'id': op_hash,
                'smirks': data['smirks'],
                'sources': ','.join(selected_sources)
            })

    # Print statistics
    print("Statistics for retained EVODEX-E operators:")
    print(f"  Number of operators: {len(retained_ops)}")
    print(f"  Sources per operator: min={min_sources}, max={max_sources}, mean={mean_sources:.2f}, median={median_sources}")

    # Save report
    report_path = os.path.join(paths['errors_dir'], 'phase3a_evodex_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write("EVODEX-E Operator Retention Report\n")
        report_file.write(f"Number of operators: {len(retained_ops)}\n")
        report_file.write(f"Sources per operator:\n")
        report_file.write(f"  min: {min_sources}\n")
        report_file.write(f"  max: {max_sources}\n")
        report_file.write(f"  mean: {mean_sources:.2f}\n")
        report_file.write(f"  median: {median_sources}\n")
        report_file.write(f"\nOperator retention breakdown:\n")
        report_file.write(f"  Operators with only 1 source (excluded): {num_ops_excluded}\n")
        report_file.write(f"  Operators retained (>=5 sources): {num_ops_retained}\n")
        report_file.write(f"  Operators with >5 sources and thus trimmed to 5: {num_ops_trimmed}\n")
        report_file.write(f"  Operators rejected for missing mapped atoms: {rejected_unmapped}\n")

    print(f"Phase 3a EVODEX-E pruning complete. Retained {len(retained_ops)} operators with >=5 sources. Up to 5 sources retained per operator.")

    # === Prune EVODEX-P to only retained P entries ===
    # Build set of surviving P IDs — up to 5 per E
    surviving_p_ids = set()
    for data in retained_ops.values():
        selected_sources = sorted(data['sources'])[:5] if data['sources'] else []
        surviving_p_ids.update(selected_sources)

    # Load EVODEX-P again
    p_rows = []
    with open(paths['evodex_p_filtered'], 'r') as pfile:
        reader = csv.DictReader(pfile)
        for row in reader:
            p_rows.append(row)

    # Filter P entries
    retained_p_rows = [row for row in p_rows if row['id'] in surviving_p_ids]

    eliminated_p_count = len(p_rows) - len(retained_p_rows)
    print(f"  EVODEX-P entries eliminated by pruning: {eliminated_p_count}")

    with open(report_path, 'a') as report_file:
        report_file.write(f"\nEVODEX-P entries eliminated by pruning: {eliminated_p_count}\n")

    # Write retained EVODEX-P
    with open(paths['evodex_p_phase3a_retained'], 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=p_rows[0].keys())
        writer.writeheader()
        writer.writerows(retained_p_rows)

    print(f"Pruned EVODEX-P: {len(retained_p_rows)} retained of {len(p_rows)} original entries.")

    # Done
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Phase 3a EVODEX-E pruning completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()