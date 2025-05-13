import csv
import os
from collections import defaultdict
from evodex.operators import extract_operator
from evodex.utils import reaction_hash
from pipeline.config import load_paths
from pipeline.version import __version__

def main():
    paths = load_paths('pipeline/config/paths.yaml')

    operator_configs = {
        'evodex_c': dict(include_stereochemistry=True, include_sigma=False, include_pi=False,
                         include_unmapped_hydrogens=True, include_unmapped_heavy_atoms=True, include_static_hydrogens=True),
        'evodex_n': dict(include_stereochemistry=True, include_sigma=True, include_pi=False,
                         include_unmapped_hydrogens=True, include_unmapped_heavy_atoms=True, include_static_hydrogens=True),
        'evodex_em': dict(include_stereochemistry=False, include_sigma=True, include_pi=True,
                          include_unmapped_hydrogens=False, include_unmapped_heavy_atoms=False, include_static_hydrogens=False),
        'evodex_cm': dict(include_stereochemistry=False, include_sigma=False, include_pi=False,
                          include_unmapped_hydrogens=False, include_unmapped_heavy_atoms=False, include_static_hydrogens=False),
        'evodex_nm': dict(include_stereochemistry=False, include_sigma=True, include_pi=False,
                          include_unmapped_hydrogens=False, include_unmapped_heavy_atoms=False, include_static_hydrogens=False),
    }

    operator_maps = {key: defaultdict(lambda: {'smirks': None, 'sources': set()}) for key in operator_configs}

    with open(paths['evodex_e'], 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            evodex_e_id = row['id']
            evodex_p = row['sources'].split(',')[0]
            smirks = row['smirks']

            for key, params in operator_configs.items():
                try:
                    op_smirks = extract_operator(smirks, **params)
                    if op_smirks.startswith(" >>") or op_smirks.endswith(">>"):
                        continue
                    op_hash = reaction_hash(op_smirks)
                    operator_maps[key][op_hash]['smirks'] = op_smirks
                    operator_maps[key][op_hash]['sources'].add(evodex_e_id)
                except Exception:
                    continue

    for key, data_map in operator_maps.items():
        out_path = paths[key]
        with open(out_path, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=['id', 'smirks', 'sources'])
            writer.writeheader()
            for rxn_hash, data in data_map.items():
                writer.writerow({'id': rxn_hash, 'smirks': data['smirks'], 'sources': ','.join(sorted(data['sources']))})

    print("Phase 4 operator mining complete.")

if __name__ == "__main__":
    main()
