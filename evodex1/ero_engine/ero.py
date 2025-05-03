

from dataclasses import dataclass
import hashlib

@dataclass(frozen=True)
class ERO:
    smirks: str
    source_ec: str  # e.g. "1.1.1.1"

    def hash(self) -> str:
        """Return a hash string uniquely identifying this ERO."""
        return hashlib.sha256(self.smirks.encode("utf-8")).hexdigest()

    def __str__(self):
        return f"ERO({self.source_ec}): {self.smirks}"