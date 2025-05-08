

import csv
import os
from collections import defaultdict
from evodex.mapping import mine_operator_from_mapping, project_operator
from evodex.utils import hash_operator
from pipeline.config import load_paths
from pipeline.version import __version__

def ensure_directories(paths: dict):
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def load_evop_data(evop_path):
    with open(evop_path, 'r') as f:
        return list(csv.DictReader(f))

def load_ero_data(ero_path):
    eros = []
    with open(ero_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            eros.append({
                'id': row['id'],
                'ero_hash': row['ero_hash'],
                'smirks': row['smirks'],
                'sources': row['sources'].split(',') if row['sources'] else []
            })
    return eros

def rank_eros_by_evop_explanation(evop_data, ero_data):
    # Map hash to ERO
    hash_to_ero = {ero['ero_hash']: ero for ero in ero_data}
    ero_explains = defaultdict(list)  # ero_hash -> list of evop idx
    for i, p in enumerate(evop_data):
        try:
            ero_obj = mine_operator_from_mapping(p['smirks'])
            ero_hash = hash_operator(ero_obj)
            if ero_hash in hash_to_ero:
                ero_explains[ero_hash].append(i)
        except Exception:
            continue
    ranked = sorted(ero_explains.items(), key=lambda x: len(x[1]), reverse=True)
    return ranked, ero_explains

def project_ero_onto_evop(evop, ero_smirks):
    # Try to project operator onto substrate/product pair
    try:
        projection = project_operator(ero_smirks, evop['substrate'], evop['product'])
        return projection['match'] if isinstance(projection, dict) and 'match' in projection else bool(projection)
    except Exception:
        return False

def assign_evop_to_best_ero(evop_data, eros):
    # For each EVODEX-P, find the best matching ERO by projection
    assignments = {}
    ero_hash_to_ero = {ero['ero_hash']: ero for ero in eros}
    for i, evop in enumerate(evop_data):
        best_ero = None
        for ero in eros:
            if project_ero_onto_evop(evop, ero['smirks']):
                best_ero = ero['ero_hash']
                break  # First matching is chosen
        if best_ero:
            assignments[evop['id']] = best_ero
        else:
            assignments[evop['id']] = None
    return assignments

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    print("Loading data...")
    evop_data = load_evop_data(paths['evodex_p'])
    ero_data = load_ero_data(paths['evodex_e_raw'])

    print("Ranking EROs by how many EVODEX-Ps they explain...")
    ranked, ero_explains = rank_eros_by_evop_explanation(evop_data, ero_data)

    print("Selecting top EROs and confirming projections...")
    # Keep EROs that explain at least one EVODEX-P
    top_eros = []
    for ero_hash, evop_idxs in ranked:
        ero = next((e for e in ero_data if e['ero_hash'] == ero_hash), None)
        if not ero:
            continue
        # Confirm at least one actual projection match (not just hash)
        any_match = False
        for idx in evop_idxs:
            if project_ero_onto_evop(evop_data[idx], ero['smirks']):
                any_match = True
                break
        if any_match:
            top_eros.append(ero)

    print(f"Assigning each EVODEX-P to the best-matching ERO...")
    evop_to_ero = assign_evop_to_best_ero(evop_data, top_eros)

    print(f"Writing {len(top_eros)} trimmed EROs to file...")
    with open(paths['evodex_e'], 'w', newline='') as out:
        fieldnames = ['id', 'ero_hash', 'smirks', 'sources']
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        for ero in top_eros:
            writer.writerow({
                'id': ero['id'],
                'ero_hash': ero['ero_hash'],
                'smirks': ero['smirks'],
                'sources': ','.join(ero['sources'])
            })

    print(f"Writing EVODEX-P to ERO mapping...")
    with open(paths['evodex_p_to_e'], 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['evop_id', 'ero_hash'])
        writer.writeheader()
        for evop_id, ero_hash in evop_to_ero.items():
            writer.writerow({'evop_id': evop_id, 'ero_hash': ero_hash or ''})

    print("Projection-based ERO trimming and mapping complete.")

if __name__ == "__main__":
    main()