from pipeline.data_preparation import main as data_preparation_main
from pipeline.operator_mining.step1_generate_partial_reactions import main as step1_main
from pipeline.operator_mining.step2_group_by_formula import main as step2_main
from pipeline.operator_mining.step3_mine_initial_eros import main as step3_main
from pipeline.operator_mining.step4_trim_and_project_eros import main as step4_main
from pipeline.web_generation.main import main as web_generation_main
from pipeline.analysis.run_analysis import main as analysis_main

def run_pipeline():
    print("Running data preparation...")
    data_preparation_main()

    print("Running mining block step 1...")
    step1_main()
    print("Running mining block step 2...")
    step2_main()
    print("Running mining block step 3...")
    step3_main()
    print("Running mining block step 4...")
    step4_main()

    print("Running website generation...")
    web_generation_main()

    print("Running EC analysis...")
    analysis_main()

    print("Pipeline execution complete.")

if __name__ == "__main__":
    run_pipeline()
