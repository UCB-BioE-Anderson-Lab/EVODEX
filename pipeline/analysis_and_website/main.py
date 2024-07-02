import os
import shutil
from pipeline.config import load_paths
from pipeline.analysis_and_website.generate_html import generate_html_pages
from pipeline.analysis_and_website.generate_svg import generate_all_svgs, generate_svgs_for_data_preparation

def setup_directories(base_dir):
    images_dir = os.path.join(base_dir, 'images')
    data_dir = os.path.join(base_dir, 'data')
    pages_dir = os.path.join(base_dir, 'pages')

    os.makedirs(images_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(pages_dir, exist_ok=True)

    return images_dir, data_dir, pages_dir

def move_csv_files(data_dir, paths):
    source_dir = os.path.join('data', 'processed')
    for filename in os.listdir(source_dir):
        if filename.endswith(".csv"):
            src_path = os.path.join(source_dir, filename)
            dest_path = os.path.join(data_dir, filename)
            shutil.copy(src_path, dest_path)
            print(f"Copied: {src_path} to {dest_path}")

    # Also copy the raw data file
    raw_src_path = paths['raw_data']
    raw_dest_path = os.path.join(data_dir, os.path.basename(raw_src_path))
    shutil.copy(raw_src_path, raw_dest_path)
    print(f"Copied: {raw_src_path} to {raw_dest_path}")

def create_website(paths):
    base_dir = 'website'
    images_dir, data_dir, pages_dir = setup_directories(base_dir)
    move_csv_files(data_dir, paths)

    ro_metadata = {
        'R': {'filename': paths['evodex_r']},
        'P': {'filename': paths['evodex_p']},
        'E': {'filename': paths['evodex_e']},
        'N': {'filename': paths['evodex_n']},
        'C': {'filename': paths['evodex_c']},
        'Em': {'filename': paths['evodex_em']},
        'Nm': {'filename': paths['evodex_nm']},
        'Cm': {'filename': paths['evodex_cm']},
        'F': {'filename': paths['evodex_f']},
        'M': {'filename': paths['evodex_m']}
    }

    data_paths = {
        'raw': paths['raw_data'],
        'filtered': paths['filtered_data'],
        'astatine': paths['astatine_data']
    }

    generate_svgs_for_data_preparation(data_paths, images_dir)

    for evodex_type, metadata in ro_metadata.items():
        generate_all_svgs(evodex_type, metadata, images_dir)

    generate_html_pages(paths, data_dir, images_dir, pages_dir, list(ro_metadata.keys()))

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    create_website(paths)
    print("Website generation complete.")

if __name__ == "__main__":
    main()
