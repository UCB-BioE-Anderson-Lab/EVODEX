# setup_directories.py
import os
import shutil

def setup_directories(base_dir, logo_path):
    os.makedirs(os.path.join(base_dir, 'images'), exist_ok=True)
    os.makedirs(os.path.join(base_dir, 'data'), exist_ok=True)
    os.makedirs(os.path.join(base_dir, 'pages'), exist_ok=True)

    # Copy the logo image to the images directory
    logo_dest = os.path.join(base_dir, 'images', 'evodex_logo.png')
    shutil.copy(logo_path, logo_dest)
    print(f"Copied logo: {logo_path} to {logo_dest}")

if __name__ == "__main__":
    base_dir = 'website'
    logo_path = 'pipeline/web_generation/templates/evodex_logo.png'
    setup_directories(base_dir, logo_path)
