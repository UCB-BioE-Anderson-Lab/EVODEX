

from typing import List
from .ero import ERO

class EROMatcher:
    def __init__(self, saved_eros: List[ERO]):
        self.saved_eros = saved_eros

    def match(self, reaction_smirks: str) -> bool:
        """
        Check if any saved ERO matches the given reaction SMIRKS.
        """
        for ero in self.saved_eros:
            if ero.smirks == reaction_smirks:
                return True
        return False