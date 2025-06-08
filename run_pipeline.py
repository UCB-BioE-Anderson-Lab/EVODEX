import subprocess
import time
import os
import requests
import gzip
import shutil
from pipeline.config import load_paths

# Define pipeline steps
PIPELINE_STEPS = [
    # "pipeline.phase1_data_preparation",
    # "pipeline.phase2_formula_pruning",
    "pipeline.phase3a_ero_mining",
    "pipeline.phase3b_ero_trimming",
    "pipeline.phase3c_ero_publishing",
    "pipeline.phase4_operator_completion",
    "pipeline.phase5_mass_subset",
    "pipeline.phase6_synthesis_subset",
    "pipeline.phase7_website",
]

def run_pipeline():
    print("=== Running full EVODEX mining pipeline ===")
    step_timings = {}

    # Step 1: Clear out data and website folders
    paths = load_paths('pipeline/config/paths.yaml')

    data_dir = 'data'
    evodex_data_dir = 'evodex/data'
    website_dir = 'website'

    # print("\n--- Clearing data and website folders ---")

    # # Clear data directory
    # os.makedirs(data_dir, exist_ok=True)
    # for filename in os.listdir(data_dir):
    #     file_path = os.path.join(data_dir, filename)
    #     if os.path.isfile(file_path) and filename.endswith('.csv'):
    #         os.remove(file_path)

    # # Clear evodex/data directory
    # os.makedirs(evodex_data_dir, exist_ok=True)
    # for filename in os.listdir(evodex_data_dir):
    #     file_path = os.path.join(evodex_data_dir, filename)
    #     if os.path.isfile(file_path) and (filename.endswith('.csv') or filename.endswith('.json')):
    #         os.remove(file_path)

    # # Clear website directory
    # os.makedirs(website_dir, exist_ok=True)
    # for root, dirs, files in os.walk(website_dir):
    #     for filename in files:
    #         file_path = os.path.join(root, filename)
    #         os.remove(file_path)

    # print("--- Data and website folders cleared ---")

    # # Step 2: Download and extract raw data
    # raw_dir = os.path.join(data_dir, 'raw')
    # os.makedirs(raw_dir, exist_ok=True)

    # print("\n--- Downloading and extracting raw data ---")
    # url = "https://github.com/hesther/enzymemap/blob/main/data/processed_reactions.csv.gz?raw=true"
    # gz_path = os.path.join(raw_dir, "processed_reactions.csv.gz")
    # csv_path = os.path.join(raw_dir, "raw_reactions.csv")

    # r = requests.get(url)
    # with open(gz_path, "wb") as f:
    #     f.write(r.content)

    # with gzip.open(gz_path, "rt") as f_in:
    #     with open(csv_path, "wt") as f_out:
    #         f_out.write(f_in.read())

    # print(f"--- Raw data downloaded and extracted to {csv_path} ---")

    for step in PIPELINE_STEPS:
        print(f"\n--- Starting {step} ---")
        start_time = time.time()
        try:
            subprocess.run(["python", "-m", step], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running {step}: {e}")
            break
        end_time = time.time()
        elapsed = end_time - start_time
        step_timings[step] = elapsed
        print(f"--- Completed {step} in {elapsed:.2f} seconds ---")

    print("\n=== Pipeline execution complete ===")
    print("\nPipeline Step Timings:")
    print("{:<50} {:>10}".format("Step", "Time (s)"))
    print("-" * 60)
    for step, elapsed in step_timings.items():
        print("{:<50} {:>10.2f}".format(step, elapsed))

    # Write timings to file
    paths = load_paths('pipeline/config/paths.yaml')
    errors_dir = paths.get('errors_dir', 'data/errors')
    os.makedirs(errors_dir, exist_ok=True)
    timings_file = os.path.join(errors_dir, "pipeline_timings.txt")

    with open(timings_file, "w") as f:
        f.write("Pipeline Step Timings:\n")
        f.write("{:<50} {:>10}\n".format("Step", "Time (s)"))
        f.write("-" * 60 + "\n")
        for step, elapsed in step_timings.items():
            f.write("{:<50} {:>10.2f}\n".format(step, elapsed))

    print(f"\nTiming report written to {timings_file}")

if __name__ == "__main__":
    run_pipeline()