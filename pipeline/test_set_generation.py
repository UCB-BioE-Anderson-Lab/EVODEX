import pandas as pd
import os
from collections import defaultdict
from evodex.utils import reaction_hash

# Paths
raw_reactions_file = os.path.join('data', 'raw', 'raw_reactions.csv')
selected_reactions_file = os.path.join('website', 'data', 'selected_reactions.csv')
output_file = os.path.join('website', 'analysis', 'test_set_reactions.csv')

# Read the input CSV files
raw_df = pd.read_csv(raw_reactions_file)
selected_df = pd.read_csv(selected_reactions_file)

# Get the set of reaction indices used in the selected reactions
selected_rxn_idx = set(selected_df['rxn_idx'])

# Initialize a dictionary to keep track of selected reaction hashes for each EC number
ec_number_hashes = defaultdict(set)

# Populate the dictionary with hashes from the selected reactions
for idx, row in selected_df.iterrows():
    ec_num = row['ec_num']
    smirks = row['mapped']
    rxn_hash = reaction_hash(smirks)
    ec_number_hashes[ec_num].add(rxn_hash)

# Initialize a dictionary to store new test reactions
new_test_reactions = defaultdict(list)

# Iterate through the raw reactions
for idx, row in raw_df.iterrows():
    ec_num = row['ec_num']
    smirks = row['mapped']
    rxn_idx = row['rxn_idx']
    
    # Skip reactions that are already in the selected set
    if rxn_idx in selected_rxn_idx:
        continue
    
    rxn_hash = reaction_hash(smirks)
    
    # Check if this reaction has a different hash and we haven't already selected two reactions for this EC number
    if rxn_hash not in ec_number_hashes[ec_num] and len(new_test_reactions[ec_num]) < 2:
        new_test_reactions[ec_num].append(row)
        ec_number_hashes[ec_num].add(rxn_hash)

# Create a DataFrame for the new test reactions
test_reactions_list = [item for sublist in new_test_reactions.values() for item in sublist]
test_reactions_df = pd.DataFrame(test_reactions_list)

# Ensure we have 2 reactions per EC number, filling with reactions from the original set if needed
for ec_num, rows in new_test_reactions.items():
    if len(rows) < 2:
        additional_reactions = raw_df[(raw_df['ec_num'] == ec_num) & (~raw_df['rxn_idx'].isin(selected_rxn_idx))]
        additional_needed = 2 - len(rows)
        for _, additional_row in additional_reactions.head(additional_needed).iterrows():
            test_reactions_df = pd.concat([test_reactions_df, pd.DataFrame([additional_row])], ignore_index=True)

# Save the test reactions to a CSV file
test_reactions_df.to_csv(output_file, index=False)

print(f"Test reactions saved to {output_file}")
