

import pandas as pd
from collections import defaultdict

def load_ec_levels(df):
    """Extract EC levels from a dataframe."""
    ec_levels = defaultdict(set)
    for ec in df['ec_num'].dropna():
        parts = ec.split('.')
        for i in range(1, 5):
            if len(parts) >= i:
                ec_levels[i].add('.'.join(parts[:i]))
    return ec_levels

def compute_coverage(raw_levels, selected_levels):
    """Compute percentage coverage for each EC level."""
    coverage = []
    for level in range(1, 5):
        raw = raw_levels[level]
        selected = selected_levels[level]
        common = raw & selected
        percent = (len(common) / len(raw)) * 100 if raw else 0
        coverage.append((level, len(selected), len(raw), len(common), percent))
    return coverage

def main():
    raw_df = pd.read_csv("data/raw/raw_reactions.csv")
    selected_df = pd.read_csv("evodex/data/selected_reactions.csv")

    raw_levels = load_ec_levels(raw_df)
    selected_levels = load_ec_levels(selected_df)

    coverage_stats = compute_coverage(raw_levels, selected_levels)

    result_df = pd.DataFrame(coverage_stats, columns=[
        'EC Level', 'Selected Count', 'Raw Count', 'Overlap Count', 'Coverage (%)'
    ])

    output_path = "data/processed/ec_coverage_report.csv"
    result_df.to_csv(output_path, index=False)
    print("Coverage report saved to:", output_path)
    print(result_df)

if __name__ == "__main__":
    main()