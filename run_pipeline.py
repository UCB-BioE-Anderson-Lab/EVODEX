from pipeline.data_preparation import main as data_preparation_main
from pipeline.mining_block import main as mining_block_main
from pipeline.analysis_and_website.main import main as analysis_and_website_main

def run_pipeline():
    print("Running data preparation...")
    data_preparation_main()

    print("Running mining block...")
    mining_block_main()

    print("Running analysis and website generation...")
    analysis_and_website_main()

    print("Pipeline execution complete.")

if __name__ == "__main__":
    run_pipeline()
