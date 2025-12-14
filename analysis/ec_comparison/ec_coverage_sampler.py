import pandas as pd
import random
from evodex.decofactor import remove_cofactors

# File path
raw_data_path = "data/raw/raw_reactions.csv"

# Load raw data
raw_df = pd.read_csv(raw_data_path)

# Pre-index EC hierarchy
ec_hierarchy = {}
for _, row in raw_df.iterrows():
    rxn_id = str(row['rxn_idx'])
    ec = row['ec_num']
    if not ec or pd.isna(ec): continue
    levels = ec.split('.')
    ec_hierarchy[rxn_id] = {
        'ec_num': ec,
        'level_1': levels[0],
        'level_2': '.'.join(levels[:2]) if len(levels) > 1 else '',
        'level_3': '.'.join(levels[:3]) if len(levels) > 2 else '',
        'level_4': '.'.join(levels) if len(levels) == 4 else ''
    }

# Initialize
num_reactions = 1000
random.seed(42)
covered_levels = {'level_1': set(), 'level_2': set(), 'level_3': set(), 'level_4': set()}
sampled_rows = []

shuffled = raw_df.sample(frac=1, random_state=42)

for _, row in shuffled.iterrows():
    if len(sampled_rows) >= num_reactions:
        break

    rxn_id = str(row['rxn_idx'])
    mapped = row['mapped']
    ec = ec_hierarchy.get(rxn_id, {}).get('ec_num')
    if not ec:
        continue

    levels = ec_hierarchy[rxn_id]
    new = False
    for lvl in ['level_1', 'level_2', 'level_3', 'level_4']:
        if levels[lvl] and levels[lvl] not in covered_levels[lvl]:
            new = True
            covered_levels[lvl].add(levels[lvl])
            break

    if not new:
        continue

    try:
        stripped = remove_cofactors(mapped)
        subs, prods = stripped.split(">>")
        if not subs or not prods:
            continue
        if '.' in subs or '.' in prods:
            continue
        from rdkit import Chem
        sub_mol = Chem.MolFromSmiles(subs)
        prod_mol = Chem.MolFromSmiles(prods)

        for atom in sub_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        for atom in prod_mol.GetAtoms():
            atom.SetAtomMapNum(0)

        sub_clean = Chem.MolToSmiles(sub_mol, isomericSmiles=True)
        prod_clean = Chem.MolToSmiles(prod_mol, isomericSmiles=True)
        smiles = f"{sub_clean}>>{prod_clean}"

        sampled_rows.append({
            'id': len(sampled_rows),
            'rxn_id': rxn_id,
            'original_rxn': mapped,
            'decofactored': stripped,
            'ec_num': ec,
            'smiles': smiles
        })
    except Exception as e:
        continue

# Save to CSV
output_path = "data/processed/ec_sampled_reactions.csv"
pd.DataFrame(sampled_rows).to_csv(output_path, index=False)
print(f"\nSaved sampled reactions to: {output_path}")