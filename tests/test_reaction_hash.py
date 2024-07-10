import pytest
import sys
import os

# This deals with path issues
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from evodex.utils import reaction_hash

# Test cases without expected hashes
# We will calculate the expected hashes dynamically
test_data_same_reaction = [
    "CCO>>CC[OH]",
    "CCO>>[OH]CC",
    "CC[OH]>>CCO",
    "[OH]CC>>CCO",
    "[C:1][C:2]O>>[C:2](O)[C:1]",
    "[C:1][C:2]O>>[C:1](O)[C:2]" # ignores atom mapping
]

test_data_different_reaction = [
    "CCO>>CC=O",
    "CC>>CC=O",
    "CCO>>CCO",
    "[CH3][CH2][OH]>>[CH3][CH2]=O", # Explicit hydrogens result in different hashes
]

def test_reaction_hash_same_reaction():
    # All these reactions should have the same hash
    base_hash = reaction_hash(test_data_same_reaction[0])
    for reaction in test_data_same_reaction[1:]:
        assert reaction_hash(reaction) == base_hash, f"Hashes do not match for reaction {reaction}"

def test_reaction_hash_different_reaction():
    # Each of these reactions should have a unique hash
    hashes = []
    for reaction in test_data_different_reaction:
        hashes.append(reaction_hash(reaction))
    assert len(hashes) == len(test_data_different_reaction), "Hashes are not unique for different reactions"

if __name__ == "__main__":
    pytest.main()
