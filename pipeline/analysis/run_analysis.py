from pipeline.analysis.aggregate_data import main as aggregate_data_main
from pipeline.analysis.generate_hierarchy import main as generate_hierarchy_main
from pipeline.analysis.analyze_hierarchy import main as analyze_hierarchy_main
from pipeline.analysis.correlation_analysis import analyze

def main():
    # Call the main method of aggregate_data.py and get the ec_map object
    ec_map = aggregate_data_main()
    
    # Call the main method of generate_hierarchy.py with ec_map
    hierarchy = generate_hierarchy_main(ec_map)
    
    # Call the main method of analyze_hierarchy.py with hierarchy
    analysis = analyze_hierarchy_main(hierarchy)
    
    # Call the main method of correlation_analysis.py with hierarchy object
    correlation_results = analyze(ec_map)
    
    return correlation_results

if __name__ == "__main__":
    correlation_results = main()
    print("Correlation analysis has been performed and results are available as an object.")
