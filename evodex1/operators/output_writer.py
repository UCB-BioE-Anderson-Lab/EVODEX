import csv
from typing import Dict, List

class OutputWriter:
    def __init__(self, output_dir: str):
        self.output_dir = output_dir

    def write_csv(self, operator_dict: Dict[str, List[str]]):
        """
        Write each list of operators (by type) to a separate CSV file in output_dir.
        Each operator is written with a header row.
        """
        for op_type, smirks_list in operator_dict.items():
            file_path = f"{self.output_dir}/EVODEX-{op_type}_reaction_operators.csv"
            with open(file_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(["smirks"])  # header
                for smirks in smirks_list:
                    writer.writerow([smirks])
