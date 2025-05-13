import csv
import os
from collections import defaultdict, Counter
from evodex.operators import extract_operator
from evodex.utils import reaction_hash
from pipeline.config import load_paths
from pipeline.version import __version__

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    print("Starting Phase 3 EVODEX-E pruning...")
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    operator_map = defaultdict(lambda: {'smirks': None, 'sources': set()})

    error_log_path = os.path.join(paths['errors_dir'], 'phase3_ero_errors.csv')
    with open(error_log_path, 'w', newline='') as errorfile:
        err_writer = csv.DictWriter(errorfile, fieldnames=['id', 'smirks', 'error_message'])
        err_writer.writeheader()

        with open(paths['evodex_p_filtered'], 'r') as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                if reader.line_num % 100 == 0:
                    print(f"Processing EVODEX-P row {reader.line_num}...")
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

    retained_ops = {k: v for k, v in operator_map.items() if len(v['sources']) >= 5}
    retained_p_ids = {pid for v in retained_ops.values() for pid in v['sources']}

    # Write evodex_e.csv
    with open(paths['evodex_e'], 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'smirks', 'sources'])
        writer.writeheader()
        for op_hash, data in retained_ops.items():
            writer.writerow({'id': op_hash, 'smirks': data['smirks'], 'sources': ','.join(sorted(data['sources']))})

    # Filter evodex_p
    with open(paths['evodex_p_filtered'], 'r') as infile, open(paths['evodex_p'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if row['id'] in retained_p_ids:
                writer.writerow(row)

    # Filter evodex_f
    retained_f_ids = set()
    with open(paths['evodex_f_filtered'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            if reader.line_num % 100 == 0:
                print(f"Checking EVODEX-F row {reader.line_num} for retention...")
            sources = row['sources'].split(',')
            if any(s in retained_p_ids for s in sources):
                retained_f_ids.add(row['id'])

    with open(paths['evodex_f_filtered'], 'r') as infile, open(paths['evodex_f'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if reader.line_num % 100 == 0:
                print(f"Writing EVODEX-F row {reader.line_num}...")
            if row['id'] in retained_f_ids:
                writer.writerow(row)

    # Write summary report
    report_lines = [
        f"EVODEX Phase 3 Final Pruning Report (version {__version__})",
        "=============================================================",
        f"Total EVODEX-E extracted: {len(operator_map)}",
        f"Retained EVODEX-E (≥5 sources): {len(retained_ops)}",
        f"Retained EVODEX-P: {len(retained_p_ids)}",
        f"Retained EVODEX-F: {len(retained_f_ids)}"
    ]

    for line in report_lines:
        print(line)

    report_path = os.path.join(paths['errors_dir'], 'Phase3_evodex_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write('\n'.join(report_lines))
    print("Phase 3 pruning complete.")

if __name__ == "__main__":
    main()