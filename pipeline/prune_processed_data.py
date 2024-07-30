import csv
import os
from collections import defaultdict

def ensure_directory(path):
    """Ensure that the directory exists."""
    if not os.path.exists(path):
        os.makedirs(path)

def read_csv(filepath):
    """Read a CSV file and return the content as a list of dictionaries."""
    with open(filepath, 'r') as infile:
        return list(csv.DictReader(infile))

def write_csv(filepath, data, fieldnames):
    """Write a list of dictionaries to a CSV file."""
    if data:  # Check if data is not empty
        with open(filepath, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)

def _parse_sources(sources):
    """
    Parse the sources field from the CSV file.

    Parameters:
    sources (str): The sources field as a string.

    Returns:
    list: A list of source strings.
    """
    sources = "" + str(sources)  # Turn the int to a string, or no effect if already a string
    sources = sources.replace('"', '')  # Remove all double quotes
    return sources.split(',')  # Split by commas

def prune_evodeex_e(data, top_n=100, max_sources=3):
    """Prune EVODEX-E data to the top N based on the count of unique EVODEX-P references, limiting sources to max_sources."""
    data_with_counts = []
    for row in data:
        sources = _parse_sources(row['sources'])
        count = len(sources)
        data_with_counts.append((row, count))

    # Sort by count in descending order and select the top N
    sorted_data = sorted(data_with_counts, key=lambda x: x[1], reverse=True)[:top_n]
    
    # Extract pruned E data and corresponding limited P references
    pruned_e = []
    pruned_p = set()
    for item in sorted_data:
        row = item[0]
        sources = _parse_sources(row['sources'])[:max_sources]  # Limit sources to max_sources
        row['sources'] = ','.join(sources)
        pruned_e.append(row)
        pruned_p.update(sources)

    print("Pruned E count:", len(pruned_e))  # Debug: print pruned E count
    print("Pruned P references:", len(pruned_p))  # Debug: print pruned P references
    
    return pruned_e, pruned_p

def prune_by_p_references(data, valid_p):
    """Prune data entries whose sources reference at least one valid EVODEX-P ID."""
    pruned_data = []
    for row in data:
        sources = _parse_sources(row['sources'])
        valid_sources = [s for s in sources if s in valid_p]
        if valid_sources:
            row['sources'] = ','.join(valid_sources)
            pruned_data.append(row)
    print(f"Pruned data count: {len(pruned_data)}")  # Debug: print pruned data count
    return pruned_data

def prune_evodeex_r(data, valid_r):
    """Prune EVODEX-R data to include only entries with valid sources and return corresponding rxn_idx."""
    pruned_r = []
    valid_rxn_idx = set()
    for row in data:
        sources = _parse_sources(row['sources'])
        if any(source in valid_r for source in sources):
            pruned_r.append(row)
            valid_rxn_idx.update(sources)
    print(f"Pruned EVODEX-R count: {len(pruned_r)}")  # Debug: print pruned EVODEX-R count
    return pruned_r, valid_rxn_idx

def main():
    input_dir = "data/processed"
    output_dir = "data/pruned_processed_data"
    ensure_directory(output_dir)

    # Read all necessary files
    evodex_e = read_csv(os.path.join(input_dir, "EVODEX-E_reaction_operators.csv"))
    evodex_p = read_csv(os.path.join(input_dir, "EVODEX-P_partial_reactions.csv"))
    evodex_r = read_csv(os.path.join(input_dir, "EVODEX-R_full_reactions.csv"))
    selected_reactions = read_csv(os.path.join(input_dir, "selected_reactions.csv"))
    
    # Additional EVODEX-* files
    additional_files = ["EVODEX-C_reaction_operators.csv", "EVODEX-Cm_reaction_operators.csv", 
                        "EVODEX-E_synthesis_subset.csv", "EVODEX-Em_reaction_operators.csv",
                        "EVODEX-F_unique_formulas.csv", "EVODEX-M_mass_spec_subset.csv", 
                        "EVODEX-M_unique_masses.csv", "EVODEX-N_reaction_operators.csv",
                        "EVODEX-Nm_reaction_operators.csv"]

    # Prune EVODEX-E
    pruned_e, valid_p = prune_evodeex_e(evodex_e, top_n=100, max_sources=3)

    # Prune EVODEX-P based on valid EVODEX-P IDs
    pruned_p = [row for row in evodex_p if row['id'] in valid_p]

    # Prune additional EVODEX-* files based on pruned EVODEX-P
    pruned_counts = {}
    for filename in additional_files:
        data = read_csv(os.path.join(input_dir, filename))
        pruned_data = prune_by_p_references(data, valid_p)
        if pruned_data:
            write_csv(os.path.join(output_dir, filename), pruned_data, pruned_data[0].keys())
            pruned_counts[filename] = len(pruned_data)

    # Prune EVODEX-R based on sources in pruned EVODEX-P
    valid_r = {source for row in pruned_p for source in _parse_sources(row['sources'])}
    pruned_r, valid_rxn_idx = prune_evodeex_r(evodex_r, valid_r)

    # Debug: Check for valid_rxn_idx and valid sources in EVODEX-R
    print(f"Valid sources in EVODEX-R: {valid_r}")
    print(f"Number of pruned EVODEX-R entries: {len(pruned_r)}")
    print(f"Valid rxn_idx: {valid_rxn_idx}")

    # Prune selected_reactions based on sources in pruned EVODEX-R
    pruned_selected_reactions = [row for row in selected_reactions if row['rxn_idx'] in valid_rxn_idx]

    # Write pruned data to new directory
    write_csv(os.path.join(output_dir, "EVODEX-E_reaction_operators.csv"), pruned_e, evodex_e[0].keys())
    write_csv(os.path.join(output_dir, "EVODEX-P_partial_reactions.csv"), pruned_p, evodex_p[0].keys())
    write_csv(os.path.join(output_dir, "EVODEX-R_full_reactions.csv"), pruned_r, evodex_r[0].keys())
    write_csv(os.path.join(output_dir, "selected_reactions.csv"), pruned_selected_reactions, selected_reactions[0].keys())

    # Print the number of entries written to each file
    print(f"EVODEX-E entries written: {len(pruned_e)}")
    print(f"EVODEX-P entries written: {len(pruned_p)}")
    print(f"EVODEX-R entries written: {len(pruned_r)}")
    print(f"selected_reactions entries written: {len(pruned_selected_reactions)}")
    for filename, count in pruned_counts.items():
        print(f"{filename} entries written: {count}")

    print(f"Pruned data written to: {output_dir}")

if __name__ == "__main__":
    main()
