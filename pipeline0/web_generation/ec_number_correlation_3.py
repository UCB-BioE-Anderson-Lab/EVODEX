import os
import pandas as pd
import json
from jinja2 import Environment, FileSystemLoader

def load_raw_data(data_dir):
    """Loads raw reactions and evodex data files from specified directory."""
    print("Loading data...")
    raw_reactions = pd.read_csv(f'{data_dir}/raw_reactions.csv')
    evodex_files = {
        'C': 'EVODEX-C_reaction_operators.csv',
        'Cm': 'EVODEX-Cm_reaction_operators.csv',
        'E': 'EVODEX-E_reaction_operators.csv',
        'Em': 'EVODEX-Em_reaction_operators.csv',
        'N': 'EVODEX-N_reaction_operators.csv',
        'Nm': 'EVODEX-Nm_reaction_operators.csv',
        'R': 'EVODEX-R_full_reactions.csv',
        'P': 'EVODEX-P_partial_reactions.csv',
        'M': 'EVODEX-M_unique_masses.csv',
        'F': 'EVODEX-F_unique_formulas.csv'
    }

    evodex_data = {}
    for evodex_type, file_name in evodex_files.items():
        evodex_data[evodex_type] = pd.read_csv(f'{data_dir}/{file_name}')
        print(f'Loaded {evodex_type} data with {evodex_data[evodex_type].shape[0]} rows')
    
    return raw_reactions, evodex_data

def expand_source_column(df, source_column='sources', delimiter=','):
    """Expands source column into separate rows for better processing."""
    df[source_column] = df[source_column].astype(str).str.replace('"', '')
    df_expanded = df.assign(**{source_column: df[source_column].str.split(delimiter)}).explode(source_column)
    return df_expanded

def split_ec_number(ec_num):
    """Splits EC number into its component parts."""
    return ec_num.split('.')

def map_evodex_to_ec(raw_reactions, evodex_data):
    ec_map = {'R': {}, 'P': {}, 'C': {}, 'Cm': {}, 'E': {}, 'Em': {}, 'N': {}, 'Nm': {}, 'M': {}, 'F': {}}
    
    # Ensure rxn_idx is treated as strings
    raw_reactions['rxn_idx'] = raw_reactions['rxn_idx'].astype(str)
    ec_dict = raw_reactions.set_index('rxn_idx')['ec_num'].to_dict()

    # Print keys and types to verify
    for key in ec_dict:
        print(key, ", ", ec_dict[key], type(key))
    
    # Process R
    evodex_r_expanded = expand_sources(evodex_data['R'], 'sources')
    for _, row in evodex_r_expanded.iterrows():
        evodex_id = str(row['id'])
        source_ids = str(row['sources']).split(',')
        for source_id in source_ids:
            source_id = source_id.strip()
            print("source_id", source_id)
            ec_num = ec_dict.get(source_id, None)
            if ec_num:
                print("found EC:", ec_num)
                if evodex_id not in ec_map['R']:
                    ec_map['R'][evodex_id] = set()
                ec_map['R'][evodex_id].add(ec_num)
            else:
                print("No EC available")
    print(f"Mapped EVODEX-R: {len(ec_map['R'])} entries")

    # Process P
    evodex_p_expanded = expand_sources(evodex_data['P'], 'sources')
    for _, row in evodex_p_expanded.iterrows():
        evodex_id = str(row['id'])
        source_ids = str(row['sources']).split(',')
        ec_set = set()
        for source_id in source_ids:
            source_id = source_id.strip()
            ec_set.update(ec_map['R'].get(source_id, set()))
        ec_map['P'][evodex_id] = ec_set
    print(f"Mapped EVODEX-P: {len(ec_map['P'])} entries")

    # Process other EVODEX types
    for evodex_type in ['C', 'Cm', 'E', 'Em', 'N', 'Nm', 'M', 'F']:
        evodex_expanded = expand_sources(evodex_data[evodex_type], 'sources')
        for _, row in evodex_expanded.iterrows():
            evodex_id = str(row['id'])
            source_ids = str(row['sources']).split(',')
            ec_set = set()
            for source_id in source_ids:
                source_id = source_id.strip()
                ec_set.update(ec_map['P'].get(source_id, set()))
            ec_map[evodex_type][evodex_id] = ec_set
        print(f"Mapped EVODEX-{evodex_type}: {len(ec_map[evodex_type])} entries")

    return ec_map

def expand_sources(df, source_column='sources', delimiter=','):
    df[source_column] = df[source_column].astype(str).str.replace('"', '')
    df_expanded = df.assign(**{source_column: df[source_column].str.split(delimiter)}).explode(source_column)
    return df_expanded

def prepare_hierarchy(ec_map, output_file):
    # Create the hierarchical dictionary structure
    hierarchy = {}

    # Function to add EC numbers and evodex IDs to the hierarchy
    def add_to_hierarchy(levels, evodex_type, evodex_id):
        if levels[0] not in hierarchy:
            hierarchy[levels[0]] = {}
        level1 = hierarchy[levels[0]]

        if levels[1] not in level1:
            level1[levels[1]] = {}
        level2 = level1[levels[1]]

        if levels[2] not in level2:
            level2[levels[2]] = {}
        level3 = level2[levels[2]]

        if levels[3] not in level3:
            level3[levels[3]] = {}
        level4 = level3[levels[3]]

        if "evodex_types" not in level4:
            level4["evodex_types"] = {}
        if evodex_type not in level4["evodex_types"]:
            level4["evodex_types"][evodex_type] = set()
        level4["evodex_types"][evodex_type].add(evodex_id)

    # Populate the hierarchical structure
    for evodex_type, evodex_dict in ec_map.items():
        print("working on: ", evodex_type)
        for evodex_id, ec_nums in evodex_dict.items():
            for ec_num in ec_nums:
                levels = split_ec_number(ec_num)
                add_to_hierarchy(levels, evodex_type, evodex_id)

    # Convert sets to lists for JSON serialization
    def convert_sets_to_lists(node):
        if isinstance(node, dict):
            if "evodex_types" in node:
                for evodex_type in node["evodex_types"]:
                    node["evodex_types"][evodex_type] = list(node["evodex_types"][evodex_type])
            for key in node:
                convert_sets_to_lists(node[key])

    convert_sets_to_lists(hierarchy)

    # Save the hierarchy to a JSON file
    with open(output_file, 'w') as f:
        json.dump(hierarchy, f, indent=4)
    print(f"Hierarchy saved to {output_file}")

    # Print summary details for assessment
    print_hierarchy_details(hierarchy)

def print_hierarchy_details(hierarchy):
    if "1" in hierarchy:
        level1 = hierarchy["1"]
        print("length of level1: ", len(level1))
    else:
        print("Key '1' not found at level 1")
        return

    if "1" in level1:
        level2 = level1["1"]
        print("length of level2: ", len(level2))
    else:
        print("Key '1' not found at level 2")
        return

    if "1" in level2:
        level3 = level2["1"]
        print("length of level3: ", len(level3))
    else:
        print("Key '1' not found at level 3")
        return

    if "1" in level3:
        level4 = level3["1"]
        print("Details at level 4: ", json.dumps(level4, indent=4))
    else:
        print("Key '1' not found at level 4")



def inject_json_to_template(json_data, template_path, output_file):
    """Injects JSON data into an HTML template using Jinja2."""
    env = Environment(loader=FileSystemLoader(os.path.dirname(template_path)))
    template = env.get_template(os.path.basename(template_path))
    rendered_html = template.render(data=json.dumps(json_data))
    with open(output_file, 'w') as output_html_file:
        output_html_file.write(rendered_html)
    print(f"HTML page saved to {output_file}")

def main():
    data_dir = 'website/data'
    output_dir = 'website/analysis'
    os.makedirs(output_dir, exist_ok=True)

    raw_reactions, evodex_data = load_raw_data(data_dir)
    ec_map = map_evodex_to_ec(raw_reactions, evodex_data)
    hierarchy_file = os.path.join(output_dir, 'ec_hierarchy.json')
    prepare_hierarchy(ec_map, hierarchy_file)
        
    with open(hierarchy_file, 'r') as f:
        json_data = json.load(f)
    
    template_path = 'pipeline/web_generation/templates/ec_hierarchy_template.html'
    output_file = os.path.join('website', 'ec_hierarchy_visualization.html')
    inject_json_to_template(json_data, template_path, output_file)

if __name__ == "__main__":
    main()
