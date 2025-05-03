import json
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import pearsonr, spearmanr
from matplotlib import rcParams
import svgutils.transform as sg
from svgutils.compose import Figure, Panel, Text

rcParams.update({'figure.autolayout': True})

def analyze(ec_map):
    correlation_analysis(ec_map)
    plot(ec_map)
    create_graphic()

def plot(ec_map):
    for evodex_type, mappings in ec_map.items():
        # Collect unique EC numbers and EVODEX IDs
        ec_numbers = set()
        evodex_ids = []
        ec_dict = {}
        for evodex_id, ec_list in mappings.items():
            evodex_ids.append(evodex_id)
            ec_dict[evodex_id] = ec_list
            for ec in ec_list:
                ec_numbers.add(ec)
        
        # Sort EC numbers
        ec_numbers = sorted(ec_numbers)
        
        # Sort EVODEX IDs alphabetically
        evodex_ids.sort()
        
        # Sort EVODEX IDs by order of EC occurrence
        evodex_ids.sort(key=lambda x: [ec_numbers.index(ec) for ec in ec_dict[x]])

        # Create a matrix for the plot
        matrix = np.zeros((len(evodex_ids), len(ec_numbers)))

        # Fill the matrix
        for i, evodex_id in enumerate(evodex_ids):
            for ec in mappings[evodex_id]:
                j = ec_numbers.index(ec)
                matrix[i, j] = 1

        # Plotting
        fig, ax = plt.subplots(figsize=(3, 9))
        ax.imshow(matrix, cmap='Greys', aspect='auto')
        
        # Set axis labels
        ax.set_xticks(np.arange(len(ec_numbers)))
        ax.set_yticks(np.arange(len(evodex_ids)))
        ax.set_xticklabels(ec_numbers, rotation=90, fontsize=1)  # Smaller font size
        y_labels = [label.replace('EVODEX-', '') for label in evodex_ids]  # Remove 'EVODEX-' prefix
        ax.set_yticklabels(y_labels, fontsize=1)  # Smaller font size
        
        # Set labels and title
        ax.set_xlabel('EC Numbers')
        ax.set_ylabel('EVODEX IDs')
        ax.set_title(f'EVODEX Type {evodex_type}')

        # Save the plot as an SVG file
        output_dir = 'website/analysis'
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(f'{output_dir}/{evodex_type}_ec_map.svg', bbox_inches='tight')
        plt.close()

def correlation_analysis(ec_map):
    results = {}
    
    for evodex_type, mappings in ec_map.items():
        # Collect unique EC numbers and EVODEX IDs
        ec_numbers = set()
        evodex_ids = []
        ec_dict = {}
        for evodex_id, ec_list in mappings.items():
            evodex_ids.append(evodex_id)
            ec_dict[evodex_id] = ec_list
            for ec in ec_list:
                ec_numbers.add(ec)
        
        # Sort EC numbers
        ec_numbers = sorted(ec_numbers)
        
        # Sort EVODEX IDs alphabetically
        evodex_ids.sort()
        
        # Sort EVODEX IDs by order of EC occurrence
        evodex_ids.sort(key=lambda x: [ec_numbers.index(ec) for ec in ec_dict[x]])

        # Prepare data for correlation calculation
        x_vals = []
        y_vals = []
        for i, evodex_id in enumerate(evodex_ids):
            for ec in ec_dict[evodex_id]:
                x_vals.append(ec_numbers.index(ec))
                y_vals.append(i)

        # Calculate Pearson and Spearman correlations
        pearson_corr, _ = pearsonr(x_vals, y_vals)
        spearman_corr, _ = spearmanr(x_vals, y_vals)
        
        # Save results
        results[evodex_type] = {
            'Pearson Correlation': pearson_corr,
            'Spearman Correlation': spearman_corr
        }

    # Save results to file
    output_dir = 'website/analysis'
    os.makedirs(output_dir, exist_ok=True)
    with open(f'{output_dir}/correlation_results.json', 'w') as f:
        json.dump(results, f, indent=4)
    
    print("Correlation analysis completed and results saved to 'website/analysis/correlation_results.json'")

# List of image files
image_files = [
    'website/analysis/R_ec_map.svg',
    'website/analysis/P_ec_map.svg',
    'website/analysis/E_ec_map.svg',
    'website/analysis/N_ec_map.svg',
    'website/analysis/C_ec_map.svg',
    'website/analysis/Em_ec_map.svg',
    'website/analysis/Nm_ec_map.svg',
    'website/analysis/Cm_ec_map.svg',
    'website/analysis/F_ec_map.svg',
    'website/analysis/M_ec_map.svg'
]

def create_graphic():
    # Layout positions for each SVG file in the combined graphic
    num_panels = len(image_files)
    panel_width = 200
    spacing = 10
    total_width = num_panels * (panel_width + spacing)
    
    # Create a list of x-coordinates with spacing between panels
    positions = [(i * (panel_width + spacing), 0) for i in range(num_panels)]

    # Load the SVG files and position them
    panels = [Panel(sg.fromfile(image_file).getroot()).move(x, y)
              for (x, y), image_file in zip(positions, image_files)]

    # Create the combined SVG figure
    fig = Figure(f"{total_width}", "900", *panels)
    output_path = 'website/analysis/combined_ec_maps.svg'
    fig.save(output_path)
    
    print(f"Combined figure saved to {output_path}")

# # Load the data
# file_path = 'website/analysis/ec_map.json'
# with open(file_path, 'r') as f:
#     ec_map = json.load(f)

# # Run the analysis
# analyze(ec_map)
