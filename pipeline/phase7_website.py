import os
import sys
import csv
from pathlib import Path

# Phase 7: Website Generation
# This phase generates the full EVODEX website. For each operator type (R, P, F, E, C, N, EM, CM, NM, M),
# it reads the corresponding CSV file, generates SVGs for each entry, and builds an HTML index page.
# All output is written to website/EVODEX_{TYPE}/ folders.

sys.path.append(os.path.join(os.path.dirname(__file__), "web_generation"))
import pipeline.web_generation.generate_svg as generate_svg
import pipeline.web_generation.generate_html as generate_html
import pipeline.web_generation.generate_css as generate_css

from pipeline.config.load_paths import load_paths

paths = load_paths('pipeline/config/paths.yaml')

# Define operator types to process
operator_types = ['R', 'P', 'F', 'E', 'C', 'N', 'Em', 'Cm', 'Nm', 'M']

# Map operator types to paths.yaml keys
operator_path_keys = {
    'R': 'evodex_r',
    'P': 'evodex_p',
    'F': 'evodex_f',
    'E': 'evodex_e',
    'C': 'evodex_c',
    'N': 'evodex_n',
    'Em': 'evodex_em',
    'Cm': 'evodex_cm',
    'Nm': 'evodex_nm',
    'M': 'evodex_m_subset',
}

# Resolve full paths using paths.yaml
operator_csv_paths = {
    op_type: paths[operator_path_keys[op_type]] for op_type in operator_types
}


# Output folder
website_root = os.path.join(os.path.dirname(__file__), '..', 'website')
os.makedirs(website_root, exist_ok=True)
generate_css.generate_css(website_root)

# Ensure images dir exists
images_dir = os.path.join(website_root, 'images')
os.makedirs(images_dir, exist_ok=True)

print("Starting Phase 7 Website Generation...")
print("=======================================")

# Process each operator type
for operator_type in operator_types:
    print(f"\nProcessing EVODEX-{operator_type}...")

    csv_path = operator_csv_paths[operator_type]
    output_dir = os.path.join(website_root, f'EVODEX_{operator_type}')

    # Read CSV
    entries = []
    with open(csv_path, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            entries.append(row)

    print(f"  Total entries: {len(entries)}")

    # Guard clause: skip SVG generation for F and M types
    if operator_type in ['F', 'M']:
        continue

    for row in entries:
        generate_svg.generate_svg(row['smirks'], f"{row['id']}.svg", images_dir)

generate_html.generate_html_pages(paths, website_root, os.path.join(website_root, 'pages'), operator_types)

print("\nPhase 7 Website Generation complete.")