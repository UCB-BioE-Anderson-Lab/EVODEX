import os
import shutil
from pipeline.config import load_paths
from pipeline.analysis_and_website.setup_directories import setup_directories
from pipeline.analysis_and_website.move_csv_files import move_and_convert_csv_files
from pipeline.analysis_and_website.generate_css import generate_css
from pipeline.analysis_and_website.generate_svg import generate_all_svgs, generate_svgs_for_data_preparation
from pipeline.analysis_and_website.generate_html import generate_html_pages
from pipeline.analysis_and_website.ec_number_correlation import main as ec_number_correlation_main

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
    move_and_convert_csv_files(data_dir, paths)
    generate_css(base_dir)

    data_paths = {
        'raw': os.path.join(data_dir, os.path.basename(paths['raw_data'])),
        'filtered': os.path.join(data_dir, os.path.basename(paths['filtered_data'])),
        'astatine': os.path.join(data_dir, os.path.basename(paths['astatine_data']))
    }

    generate_svgs_for_data_preparation(data_paths, images_dir)

    ro_metadata = {
        'R': os.path.join(data_dir, os.path.basename(paths['evodex_r'])),
        'P': os.path.join(data_dir, os.path.basename(paths['evodex_p'])),
        'E': os.path.join(data_dir, os.path.basename(paths['evodex_e'])),
        'Em': os.path.join(data_dir, os.path.basename(paths['evodex_em'])),
        'N': os.path.join(data_dir, os.path.basename(paths['evodex_n'])),
        'Nm': os.path.join(data_dir, os.path.basename(paths['evodex_nm'])),
        'C': os.path.join(data_dir, os.path.basename(paths['evodex_c'])),
        'Cm': os.path.join(data_dir, os.path.basename(paths['evodex_cm']))
    }

    for evodex_type, csv_path in ro_metadata.items():
        generate_all_svgs(evodex_type, csv_path, images_dir)

    generate_html_pages(paths, data_dir, images_dir, pages_dir, list(ro_metadata.keys()) + ['F', 'M'])
    ec_number_correlation_main()
    
def main():
    paths = load_paths('pipeline/config/paths.yaml')
    create_website(paths)
    print("Website generation complete.")

if __name__ == "__main__":
    main()
