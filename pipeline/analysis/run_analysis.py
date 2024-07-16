import aggregate_data
import generate_hierarchy
import analyze_hierarchy
from correlation_analysis import analyze

def main():
    # Call the main method of aggregate_data.py and get the ec_map object
    ec_map = aggregate_data.main()
    
    # Call the main method of generate_hierarchy.py with ec_map
    hierarchy = generate_hierarchy.main(ec_map)
    
    # Call the main method of analyze_hierarchy.py with hierarchy
    analysis = analyze_hierarchy.main(hierarchy)
    
    # Call the main method of correlation_analysis.py with hierarchy object
    correlation_results = analyze(ec_map)
    
    return correlation_results

if __name__ == "__main__":
    correlation_results = main()
    print("Correlation analysis has been performed and results are available as an object.")
