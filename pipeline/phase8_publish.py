import os
import shutil
from pipeline.config import load_paths

def publish_csv_files():
    # Load paths
    paths = load_paths('pipeline/config/paths.yaml')
    
    source_dir = os.path.join('data', 'processed')
    website_dir = 'website/data'
    evodex_dir = 'evodex/data'
    
    os.makedirs(website_dir, exist_ok=True)
    os.makedirs(evodex_dir, exist_ok=True)
    
    # Files to copy
    files_to_copy = [f for f in os.listdir(source_dir) if f.startswith('EVODEX-') and f.endswith('.csv')]
    files_to_copy += ['selected_reactions.csv', 'selected_data.csv']
    
    # Copy to website/data
    for filename in files_to_copy:
        src_path = os.path.join(source_dir, filename)
        if os.path.exists(src_path):
            dest_path = os.path.join(website_dir, filename)
            shutil.copy(src_path, dest_path)
            print(f"Copied to website: {src_path} -> {dest_path}")
    
    # Copy to evodex/data
    for filename in files_to_copy:
        src_path = os.path.join(source_dir, filename)
        if os.path.exists(src_path):
            dest_path = os.path.join(evodex_dir, filename)
            shutil.copy(src_path, dest_path)
            print(f"Copied to evodex: {src_path} -> {dest_path}")

if __name__ == "__main__":
    print("Starting Phase 8: Publishing EVODEX data...")
    publish_csv_files()
    print("Phase 8 complete: Data published to website/data and evodex/data.")