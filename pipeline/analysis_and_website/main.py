# main.py
import os
import shutil
from pipeline.config import load_paths
from pipeline.analysis_and_website.setup_directories import setup_directories
from pipeline.analysis_and_website.move_csv_files import move_csv_files
from pipeline.analysis_and_website.generate_css import generate_css
from pipeline.analysis_and_website.generate_svg import generate_all_svgs, generate_svgs_for_data_preparation
from pipeline.analysis_and_website.generate_html import generate_html_pages

def clean_website_directory(base_dir):
    if os.path.exists(base_dir):
        shutil.rmtree(base_dir)
        print(f"Deleted all contents in {base_dir}")
    os.makedirs(base_dir)
    print(f"Created empty directory {base_dir}")

def create_website(paths):
    base_dir = 'website'
    clean_website_directory(base_dir)
    logo_path = paths['logo_image']
    images_dir, data_dir, pages_dir = setup_directories(base_dir, logo_path)
    move_csv_files(data_dir, paths)
    generate_css(base_dir)

    ro_metadata = {
        'R': {'filename': paths['evodex_r']},
        'P': {'filename': paths['evodex_p']},
        'E': {'filename': paths['evodex_e']},
        'Em': {'filename': paths['evodex_em']},
        'N': {'filename': paths['evodex_n']},
        'Nm': {'filename': paths['evodex_nm']},
        'C': {'filename': paths['evodex_c']},
        'Cm': {'filename': paths['evodex_cm']}
    }

    data_paths = {
        'raw': paths['raw_data'],
        'filtered': paths['filtered_data'],
        'astatine': paths['astatine_data']
    }

    generate_svgs_for_data_preparation(data_paths, images_dir)

    for evodex_type, metadata in ro_metadata.items():
        generate_all_svgs(evodex_type, metadata, images_dir)

    generate_html_pages(paths, data_dir, images_dir, pages_dir, list(ro_metadata.keys()) + ['F', 'M'])

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    create_website(paths)
    print("Website generation complete.")

if __name__ == "__main__":
    main()
