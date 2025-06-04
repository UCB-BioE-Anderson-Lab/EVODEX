from pipeline.phase1_data_preparation import main as data_preparation_main
from pipeline.phase2_mining_block import main as mining_block_main
from pipeline.web_generation.main import main as web_generation_main
from pipeline.analysis.run_analysis import main as analysis_main

def run_pipeline():
    print("Running data preparation...")
    data_preparation_main()

    print("Running mining block...")
    mining_block_main()

    print("Running website generation...")
    web_generation_main()

    print("Running EC analysis...")
    analysis_main()

    print("Pipeline execution complete.")

if __name__ == "__main__":
    run_pipeline()
