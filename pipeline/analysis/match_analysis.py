import pandas as pd
import os

from evodex.evaluation import match_operators
from evodex.decofactor import remove_cofactors

# Path to the input CSV file
input_file = os.path.join('website', 'analysis', 'test_set_reactions.csv')
# Path to the output TSV file
output_file = os.path.join('website', 'analysis', 'output_matches.tsv')

# Read the CSV file
df = pd.read_csv(input_file)

# Initialize lists to store results
matches = []
errors = []

# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    rxn_idx = row['rxn_idx']
    smirks = row['mapped']
    
    try:
        # Run the match_operator function
        matched_operators = match_operators(smirks, 'E')
        # Append the result to the matches list
        matches.append({'rxn_idx': rxn_idx, 'matches': matched_operators, 'errors': ''})
    except Exception as e:
        # If an error occurs, append the error message to the errors list
        errors.append({'rxn_idx': rxn_idx, 'matches': '', 'errors': str(e)})

# Combine matches and errors into a single DataFrame
results_df = pd.DataFrame(matches + errors)

# Save the results to a TSV file
results_df.to_csv(output_file, sep='\t', index=False)
