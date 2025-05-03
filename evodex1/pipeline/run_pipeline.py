import os
import shutil
from evodex1.pipeline import ero_mining
from evodex1.pipeline import data_preparation

def run_pipeline():
    input_csv = "data/processed/EVODEX-R_full_reactions.csv"
    output_dir = "evodex1_outputs"
    threshold = 3

    # Clean up and recreate output directory
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    print("Running data preparation...")
    data_preparation.main()

    print("Starting ERO mining pipeline...")
    ero_mining.main(input_csv=input_csv, output_dir=output_dir, threshold=threshold)
    print("ERO mining complete.")

    # Placeholder: call test_set_generation.main() if implemented
    # print("Generating test set...")
    # test_set_generation.main()

if __name__ == "__main__":
    run_pipeline()
