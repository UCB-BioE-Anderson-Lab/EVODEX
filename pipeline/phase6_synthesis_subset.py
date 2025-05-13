


import csv
import os
from collections import defaultdict
from evodex.decofactor import remove_cofactors
from evodex.utils import reaction_hash
from pipeline.config import load_paths
from pipeline.version import __version__

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    print("Starting Phase 6: Synthesis subset generation...")

    evodex_p_map = {}
    evodex_e_map = {}

    with open(paths['evodex_p'], 'r') as p_file:
        p_reader = csv.DictReader(p_file)
        for row in p_reader:
            rxn_hash = reaction_hash(row['smirks'])
            evodex_p_map[rxn_hash] = row['id']

    with open(paths['evodex_e'], 'r') as e_file:
        e_reader = csv.DictReader(e_file)
        for row in e_reader:
            sources = row['sources'].split(',')
            for source in sources:
                evodex_e_map.setdefault(source, []).append(row['id'])

    evodex_e_subset = set()
    evodex_e_sources = {}

    with open(paths['evodex_r'], 'r') as infile, open(paths['evodex_e_synthesis'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=['id', 'sources'])
        writer.writeheader()

        for row in reader:
            try:
                smirks = row['smirks']
                partial = remove_cofactors(smirks)
                if partial.startswith(">>") or partial.endswith(">>"):
                    continue
                rxn_hash = reaction_hash(partial)
                evodex_p_id = evodex_p_map.get(rxn_hash)
                if evodex_p_id:
                    evodex_e_ids = evodex_e_map.get(evodex_p_id)
                    if evodex_e_ids:
                        for eid in evodex_e_ids:
                            evodex_e_subset.add(eid)
                            evodex_e_sources.setdefault(eid, set()).add(evodex_p_id)
            except Exception:
                continue

        for eid in sorted(evodex_e_subset):
            writer.writerow({'id': eid, 'sources': ','.join(sorted(evodex_e_sources[eid]))})

    print("Phase 6 complete: Synthesis subset written.")

if __name__ == "__main__":
    main()