

from typing import Optional
from .ero import ERO
from evodex1.operators.operators import extract_operator

class EROGenerator:
    def __init__(self, include_stereochemistry=True):
        self.include_stereochemistry = include_stereochemistry

    def generate(self, reaction_smirks: str, source_ec: str) -> Optional[ERO]:
        """
        Generate an ERO from a reaction SMIRKS string.
        Returns None if operator extraction fails.
        """
        try:
            cleaned_smirks = extract_operator(
                reaction_smirks,
                include_stereochemistry=self.include_stereochemistry
            )
            return ERO(smirks=cleaned_smirks, source_ec=source_ec)
        except Exception as e:
            print(f"[EROGenerator] Failed to extract operator: {e}")
            return None