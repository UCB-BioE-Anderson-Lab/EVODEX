import csv
import os
from collections import defaultdict, Counter
from evodex.operators import extract_operator
from evodex.utils import reaction_hash
from pipeline.config import load_paths
from pipeline.version import __version__

# Phase 3: EVODEX-E Calculation and Final Pruning
#
# This phase extracts EVODEX-E operators from the filtered EVODEX-P set, then prunes the EVODEX-P, EVODEX-R, and EVODEX-F 
# layers to retain only those linked to well-supported EVODEX-E (≥5 sources). Final sequential EVODEX.1-* IDs are assigned 
# and the finalized data stack is written.

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    print("Starting Phase 3 EVODEX-E pruning and final ID assignment...")
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    # Load EVODEX-R preliminary entries
    r_entries = []
    with open(paths['evodex_r_preliminary'], 'r') as rfile:
        reader = csv.DictReader(rfile)
        for row in reader:
            r_entries.append(row)
    total_r = len(r_entries)

    # Assign EVODEX.1-R# IDs to R entries (sequential)
    r_id_map = {}
    for idx, entry in enumerate(r_entries, start=1):
        evodex_r_id = f"EVODEX.1-R{idx}"
        r_id_map[entry['id']] = evodex_r_id
        entry['id'] = evodex_r_id

    # Load EVODEX-P filtered entries and extract operators within error log context
    p_entries = []
    error_log_path = os.path.join(paths['errors_dir'], 'phase3_ero_errors.csv')
    with open(error_log_path, 'w', newline='') as errorfile:
        err_writer = csv.DictWriter(errorfile, fieldnames=['id', 'smirks', 'error_message'])
        err_writer.writeheader()

        with open(paths['evodex_p_filtered'], 'r') as pfile:
            reader = csv.DictReader(pfile)
            for row in reader:
                p_entries.append(row)

        total_p = len(p_entries)

        # Extract operators from EVODEX-P and build operator_map
        operator_map = defaultdict(lambda: {'smirks': None, 'sources': set()})
        for i, row in enumerate(p_entries, start=1):
            if i % 100 == 0:
                print(f"Processing EVODEX-P row {i} for operator extraction...")
            try:
                op_smirks = extract_operator(
                    row['smirks'],
                    include_stereochemistry=True,
                    include_sigma=True,
                    include_pi=True,
                    include_unmapped_hydrogens=True,
                    include_unmapped_heavy_atoms=True,
                    include_static_hydrogens=True
                )
                if op_smirks.startswith(">>") or op_smirks.endswith(">>"):
                    continue
                op_hash = reaction_hash(op_smirks)
                operator_map[op_hash]['smirks'] = op_smirks
                operator_map[op_hash]['sources'].add(row['id'])
            except Exception as e:
                err_writer.writerow({'id': row['id'], 'smirks': row['smirks'], 'error_message': str(e)})

    # Retain operators supported by >=5 EVODEX-P sources
    retained_ops = {k: v for k, v in operator_map.items() if len(v['sources']) >= 5}
    retained_p_ids_set = {pid for v in retained_ops.values() for pid in v['sources']}

    # Filter EVODEX-P entries to only those in retained_p_ids_set
    filtered_p_entries = [row for row in p_entries if row['id'] in retained_p_ids_set]

    # Assign EVODEX.1-P# IDs to retained EVODEX-P entries (sequential, sorted by descending source count)
    p_source_counts = Counter()
    for op_hash, data in retained_ops.items():
        for pid in data['sources']:
            p_source_counts[pid] += 1
    # Sort filtered_p_entries by source count descending, then by original id for stability
    filtered_p_entries.sort(key=lambda x: (-p_source_counts[x['id']], x['id']))
    p_id_map = {}
    for idx, entry in enumerate(filtered_p_entries, start=1):
        evodex_p_id = f"EVODEX.1-P{idx}"
        p_id_map[entry['id']] = evodex_p_id

    # Update 'sources' field in filtered_p_entries to reference EVODEX-R IDs
    # Original sources field contains EVODEX-R IDs (preliminary), map them to EVODEX.1-R#
    for entry in filtered_p_entries:
        original_sources = entry.get('sources', '')
        if original_sources.strip():
            source_ids = [s.strip() for s in original_sources.split(',')]
            mapped_sources = [r_id_map.get(s, s) for s in source_ids]
            entry['sources'] = ','.join(mapped_sources)
        else:
            entry['sources'] = ''

    # Load EVODEX-F filtered entries
    f_entries = []
    with open(paths['evodex_f_filtered'], 'r') as ffile:
        reader = csv.DictReader(ffile)
        for row in reader:
            f_entries.append(row)
    total_f = len(f_entries)

    # Retain EVODEX-F entries whose sources include any retained EVODEX-P ids
    retained_f_entries = []
    for i, row in enumerate(f_entries, start=1):
        if i % 100 == 0:
            print(f"Checking EVODEX-F row {i} for retention...")
        sources = row.get('sources', '')
        source_ids = [s.strip() for s in sources.split(',')] if sources.strip() else []
        if any(s in retained_p_ids_set for s in source_ids):
            retained_f_entries.append(row)
    retained_f_ids_set = {row['id'] for row in retained_f_entries}

    # Assign EVODEX.1-F# IDs to retained EVODEX-F entries (sequential, sorted by descending source count)
    f_source_counts = Counter()
    for row in retained_f_entries:
        srcs = row.get('sources', '')
        src_list = [s.strip() for s in srcs.split(',')] if srcs.strip() else []
        f_source_counts[row['id']] = len(src_list)
    retained_f_entries.sort(key=lambda x: (-f_source_counts[x['id']], x['id']))
    f_id_map = {}
    for idx, entry in enumerate(retained_f_entries, start=1):
        evodex_f_id = f"EVODEX.1-F{idx}"
        f_id_map[entry['id']] = evodex_f_id

    # Update 'sources' field in retained_f_entries to reference EVODEX-P IDs
    for entry in retained_f_entries:
        original_sources = entry.get('sources', '')
        if original_sources.strip():
            source_ids = [s.strip() for s in original_sources.split(',')]
            mapped_sources = [p_id_map.get(s, s) for s in source_ids if p_id_map.get(s, s).startswith('EVODEX.1-P')]
            entry['sources'] = ','.join(mapped_sources)
        else:
            entry['sources'] = ''

    # For EVODEX-E, keep hash IDs as IDs
    # Update 'sources' field to reference EVODEX-P IDs
    retained_ops_list = list(retained_ops.items())
    # Sort retained_ops by descending number of sources
    retained_ops_list.sort(key=lambda x: -len(x[1]['sources']))

    # Update sources field for EVODEX-E to mapped EVODEX-P IDs
    for op_hash, data in retained_ops_list:
        mapped_sources = [p_id_map.get(pid, pid) for pid in data['sources']]
        data['sources'] = ','.join(sorted(mapped_sources))
    total_e = len(operator_map)
    retained_e = len(retained_ops)

    e_id_map = {}
    for idx, (op_hash, data) in enumerate(retained_ops_list, start=1):
        evodex_e_id = f"EVODEX.1-E{idx}"
        e_id_map[op_hash] = evodex_e_id

    # Write EVODEX-R final file
    with open(paths['evodex_r'], 'w', newline='') as outfile:
        fieldnames = r_entries[0].keys() if r_entries else ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for entry in r_entries:
            writer.writerow(entry)

    # Write EVODEX-P final file
    with open(paths['evodex_p'], 'w', newline='') as outfile:
        fieldnames = filtered_p_entries[0].keys() if filtered_p_entries else ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for entry in filtered_p_entries:
            # Update id to assigned EVODEX-P ID
            original_id = entry['id']
            entry['id'] = p_id_map.get(original_id, original_id)
            writer.writerow(entry)

    # Write EVODEX-F final file
    with open(paths['evodex_f'], 'w', newline='') as outfile:
        fieldnames = retained_f_entries[0].keys() if retained_f_entries else ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for entry in retained_f_entries:
            original_id = entry['id']
            entry['id'] = f_id_map.get(original_id, original_id)
            writer.writerow(entry)

    # Write EVODEX-E final file
    with open(paths['evodex_e'], 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for op_hash, data in retained_ops_list:
            writer.writerow({'id': e_id_map.get(op_hash, op_hash), 'smirks': data['smirks'], 'sources': data['sources']})

    import statistics

    e_counts = [len(data['sources'].split(',')) if data['sources'] else 0 for _, data in retained_ops_list]

    report_lines = [
        f"EVODEX Phase 3 Final Pruning and ID Assignment Report (version {__version__})",
        "===========================================================================",
        f"Total EVODEX-R preliminary entries: {total_r}",
        f"Total EVODEX-P filtered entries (input to Phase 3): {total_p}",
        f"Total EVODEX-F filtered entries (input to Phase 3): {total_f}",
        f"Total EVODEX-E extracted (from EVODEX-P): {total_e}",
        "",
        f"Final retained EVODEX-R entries: {total_r}",
        f"Final retained EVODEX-P entries (≥5 sources): {len(filtered_p_entries)}",
        f"Final retained EVODEX-F entries (linked to retained P): {len(retained_f_entries)}",
        f"Final retained EVODEX-E entries (≥5 P sources): {retained_e}",
        "",
        f"EVODEX-E sources per operator statistics:",
        f"  Min sources: {min(e_counts) if e_counts else 0}",
        f"  Max sources: {max(e_counts) if e_counts else 0}",
        f"  Mean sources: {statistics.mean(e_counts) if e_counts else 0:.2f}",
        f"  Median sources: {statistics.median(e_counts) if e_counts else 0}",
        "",
        "ID assignment scheme:",
        "  EVODEX-R: EVODEX.1-R# (sequential)",
        "  EVODEX-P: EVODEX.1-P# (sequential, sorted by source count)",
        "  EVODEX-F: EVODEX.1-F# (sequential, sorted by source count)",
        "  EVODEX-E: EVODEX.1-E# (sequential, sorted by source count)",
    ]

    for line in report_lines:
        print(line)

    report_path = os.path.join(paths['errors_dir'], 'Phase3_evodex_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write('\n'.join(report_lines))

    print("Phase 3 pruning and final ID assignment complete.")

if __name__ == "__main__":
    main()