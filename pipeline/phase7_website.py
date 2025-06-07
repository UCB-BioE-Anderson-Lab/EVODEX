import os
import sys
import csv
from pathlib import Path
import time

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
start_time = time.time()
print("Phase 7 website generation started...")

# Define operator types to process
operator_types = ['R', 'P', 'F', 'E', 'C', 'N', 'Em', 'Cm', 'Nm', 'M']

# Map operator types to paths.yaml keys
operator_path_keys = {
    'R': 'evodex_r_published',
    'P': 'evodex_p_published',
    'F': 'evodex_f_published',
    'E': 'evodex_e_published',
    'C': 'evodex_c_published',
    'N': 'evodex_n_published',
    'Em': 'evodex_em_published',
    'Cm': 'evodex_cm_published',
    'Nm': 'evodex_nm_published',
    'M': 'evodex_m_subset_published',
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

# Copy logo image to images dir
logo_src = paths['logo_image']
logo_dest = os.path.join(images_dir, 'evodex_logo.png')
from shutil import copyfile
copyfile(logo_src, logo_dest)
print(f"Copied logo: {logo_src} -> {logo_dest}")

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
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Phase 7 website generation completed in {elapsed_time:.2f} seconds.")