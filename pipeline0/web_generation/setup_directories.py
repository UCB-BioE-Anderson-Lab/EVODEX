# setup_directories.py
import os
import shutil

def setup_directories(base_dir, logo_path):
    images_dir = os.path.join(base_dir, 'images')
    data_dir = os.path.join(base_dir, 'data')
    pages_dir = os.path.join(base_dir, 'pages')

    os.makedirs(images_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(pages_dir, exist_ok=True)

    # Copy the logo image to the images directory
    logo_dest = os.path.join(images_dir, 'evodex_logo.png')
    shutil.copy(logo_path, logo_dest)
    print(f"Copied logo: {logo_path} to {logo_dest}")

    return images_dir, data_dir, pages_dir

if __name__ == "__main__":
    base_dir = 'website'
    logo_path = 'pipeline/web_generation/templates/evodex_logo.png'
    setup_directories(base_dir, logo_path)
