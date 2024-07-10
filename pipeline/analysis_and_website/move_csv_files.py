# move_csv_files.py
import os
import shutil

def move_csv_files(data_dir, paths):
    source_dir = os.path.join('data', 'processed')
    for filename in os.listdir(source_dir):
        if filename.endswith(".csv"):
            src_path = os.path.join(source_dir, filename)
            dest_path = os.path.join(data_dir, filename)
            shutil.copy(src_path, dest_path)
            print(f"Copied: {src_path} to {dest_path}")

    raw_src_path = paths['raw_data']
    raw_dest_path = os.path.join(data_dir, os.path.basename(raw_src_path))
    shutil.copy(raw_src_path, raw_dest_path)
    print(f"Copied: {raw_src_path} to {raw_dest_path}")

if __name__ == "__main__":
    data_dir = 'website/data'
    paths = load_paths('pipeline/config/paths.yaml')
    move_csv_files(data_dir, paths)
