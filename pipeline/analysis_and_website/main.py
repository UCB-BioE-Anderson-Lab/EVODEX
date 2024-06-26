import os
import shutil
from pipeline.config import load_paths
from pipeline.analysis_and_website.generate_html import generate_html_pages
from pipeline.analysis_and_website.generate_svg import generate_all_svgs

def setup_directories(base_dir):
    images_dir = os.path.join(base_dir, 'images')
    data_dir = os.path.join(base_dir, 'data')
    pages_dir = os.path.join(base_dir, 'pages')

    os.makedirs(images_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(pages_dir, exist_ok=True)

    return images_dir, data_dir, pages_dir

def move_csv_files(data_dir, paths):
    for key, path in paths.items():
        if 'processed' in path:
            csv_file = path.split('/')[-1]
            if os.path.exists(path):
                shutil.move(path, os.path.join(data_dir, csv_file))
                # Update the path in the paths dictionary to reflect the new location
                paths[key] = os.path.join(data_dir, csv_file)
                print(f"Moved: {path} to {os.path.join(data_dir, csv_file)}")
            else:
                print(f"File not found: {path}")

def create_website(paths):
    base_dir = 'website'
    images_dir, data_dir, pages_dir = setup_directories(base_dir)
    move_csv_files(data_dir, paths)

    ro_metadata = {
        'E': {
            'filename': paths['evodex_e'],
            'title': 'Electronic Reaction Operator'
        },
        'C': {
            'filename': paths['evodex_c'],
            'title': 'Core Reaction Operator'
        },
        'N': {
            'filename': paths['evodex_n'],
            'title': 'Nearest Neighbor Reaction Operator'
        },
        'Em': {
            'filename': paths['evodex_em'],
            'title': 'Minimal Electronic Reaction Operator'
        },
        'Cm': {
            'filename': paths['evodex_cm'],
            'title': 'Minimal Core Reaction Operator'
        },
        'Nm': {
            'filename': paths['evodex_nm'],
            'title': 'Minimal Nearest Neighbor Reaction Operator'
        }
    }

    generate_all_svgs(ro_metadata, data_dir, images_dir)
    generate_html_pages(paths, data_dir, images_dir, pages_dir)

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    create_website(paths)
    print("Website generation complete.")

if __name__ == "__main__":
    main()
