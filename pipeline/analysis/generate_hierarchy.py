import json

def split_ec_number(ec_num):
    """Splits EC number into its component parts."""
    return ec_num.split('.')

def prepare_hierarchy(ec_map, output_file):
    # Create the hierarchical dictionary structure
    hierarchy = {}

    # Function to add EC numbers and evodex IDs to the hierarchy
    def add_to_hierarchy(levels, evodex_type, evodex_id):
        if levels[0] not in hierarchy:
            hierarchy[levels[0]] = {}
        level1 = hierarchy[levels[0]]

        if len(levels) > 1 and levels[1] not in level1:
            level1[levels[1]] = {}
        level2 = level1.get(levels[1], {})

        if len(levels) > 2 and levels[2] not in level2:
            level2[levels[2]] = {}
        level3 = level2.get(levels[2], {})

        if len(levels) > 3 and levels[3] not in level3:
            level3[levels[3]] = {}
        level4 = level3.get(levels[3], {})

        if "evodex_types" not in level4:
            level4["evodex_types"] = {}
        if evodex_type not in level4["evodex_types"]:
            level4["evodex_types"][evodex_type] = set()
        level4["evodex_types"][evodex_type].add(evodex_id)

        if len(levels) > 1:
            level1[levels[1]] = level2
        if len(levels) > 2:
            level2[levels[2]] = level3
        if len(levels) > 3:
            level3[levels[3]] = level4

    # Populate the hierarchical structure
    for evodex_type, evodex_dict in ec_map.items():
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

    return hierarchy

def main(ec_map):
    output_file = 'website/analysis/ec_hierarchy.json'
    hierarchy = prepare_hierarchy(ec_map, output_file)
    return hierarchy

if __name__ == "__main__":
    import aggregate_data
    ec_map = aggregate_data.main()
    hierarchy = main(ec_map)
    print("Hierarchy generated and is available as an object.")
