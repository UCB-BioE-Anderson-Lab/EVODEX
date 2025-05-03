

from collections import defaultdict
from typing import Dict, List, Set
from .ero import ERO

class ERORegistry:
    def __init__(self, promotion_threshold: int = 3):
        self.usage_count: Dict[str, int] = defaultdict(int)
        self.sources: Dict[str, List[str]] = defaultdict(list)
        self.saved_eros: Set[str] = set()
        self.promotion_threshold = promotion_threshold

    def register(self, ero: ERO):
        h = ero.hash()
        self.usage_count[h] += 1
        self.sources[h].append(ero.source_ec)
        self.promote_if_threshold(ero)

    def increment_usage(self, ero: ERO):
        h = ero.hash()
        self.usage_count[h] += 1

    def get_sources(self, ero: ERO) -> List[str]:
        return self.sources[ero.hash()]

    def is_saved(self, ero: ERO) -> bool:
        return ero.hash() in self.saved_eros

    def promote_if_threshold(self, ero: ERO) -> bool:
        h = ero.hash()
        if self.usage_count[h] >= self.promotion_threshold:
            self.saved_eros.add(h)
            return True
        return False

    def get_saved_eros(self) -> Set[str]:
        return self.saved_eros