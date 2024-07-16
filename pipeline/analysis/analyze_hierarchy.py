import json

def analyze_hierarchy(hierarchy):
    """Analyze the hierarchy to count evodex IDs by type at each node of the EC tree."""
    analysis = {}

    def count_evodex_types(node):
        """Recursively count evodex types in the hierarchy."""
        evodex_counts = {}
        
        # Base case: if there are evodex_types in the node, count them
        if 'evodex_types' in node:
            for evodex_type, ids in node['evodex_types'].items():
                if evodex_type not in evodex_counts:
                    evodex_counts[evodex_type] = 0
                evodex_counts[evodex_type] += len(ids)

        # Recursive case: count evodex types in child nodes
        if isinstance(node, dict):
            for key, child in node.items():
                if key != 'evodex_types':
                    child_counts = count_evodex_types(child)
                    for evodex_type, count in child_counts.items():
                        if evodex_type not in evodex_counts:
                            evodex_counts[evodex_type] = 0
                        evodex_counts[evodex_type] += count
        
        return evodex_counts

    def analyze_node(node, path=''):
        """Recursively analyze each node and add evodex type counts to the analysis."""
        if path not in analysis:
            analysis[path] = {}

        # Count evodex types at the current node
        evodex_counts = count_evodex_types(node)
        analysis[path] = evodex_counts

        # Analyze child nodes
        if isinstance(node, dict):
            for key, child in node.items():
                if key != 'evodex_types':
                    new_path = f"{path}.{key}" if path else key
                    analyze_node(child, new_path)

    # Start the analysis from the root of the hierarchy
    analyze_node(hierarchy)

    return analysis

def main(hierarchy):
    analysis = analyze_hierarchy(hierarchy)
    
    output_file = 'website/analysis/ec_hierarchy_analysis.json'
    with open(output_file, 'w') as f:
        json.dump(analysis, f, indent=4)
    print(f"Analysis saved to {output_file}")
    
    return analysis

if __name__ == "__main__":
    import generate_hierarchy
    ec_map = generate_hierarchy.main()
    hierarchy = generate_hierarchy.prepare_hierarchy(ec_map, 'website/analysis/ec_hierarchy.json')
    analysis = main(hierarchy)
    print("Analysis generated and is available as an object.")
