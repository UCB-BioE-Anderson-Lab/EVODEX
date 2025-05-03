

from evodex1.ero_engine.registry import ERORegistry
from evodex1.ero_engine.generator import EROGenerator
from evodex1.ero_engine.matcher import EROMatcher
from evodex1.ero_engine.ero import ERO
from evodex1.operators.synthesizer import OperatorSynthesizer
from evodex1.operators.output_writer import OutputWriter
import os
import shutil
import csv
from typing import List

def main(input_csv: str, output_dir: str, threshold: int = 3):
    # Ensure the output directory is fresh
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Initialize registry and generator
    registry = ERORegistry(promotion_threshold=threshold)
    generator = EROGenerator()
    all_eros = []

    # Read reactions from input CSV and process
    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ec = row.get("ec_num", "unknown")
            smirks = row.get("smirks")
            if not smirks:
                continue

            # Step 1: Try to match saved EROs
            saved_smirks = registry.get_saved_eros()
            matcher = EROMatcher([ERO(smirks=s, source_ec="") for s in saved_smirks])
            if matcher.match(smirks):
                continue  # already explained by a saved ERO

            # Step 2: Attempt new ERO generation
            ero = generator.generate(smirks, source_ec=ec)
            if ero:
                registry.register(ero)
                all_eros.append(ero)

    # Step 3: Aggregate higher-level operators
    synthesizer = OperatorSynthesizer(all_eros)
    operator_dict = synthesizer.from_eros()

    # Step 4: Output to files
    writer = OutputWriter(output_dir)
    writer.write_csv(operator_dict)