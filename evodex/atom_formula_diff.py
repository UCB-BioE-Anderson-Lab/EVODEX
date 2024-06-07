import csv
import hashlib
from typing import Dict, Optional
from rdkit import Chem
from rdkit.Chem import AllChem

def parse_smirks(smirks: str) -> Dict[str, int]:
    """
    Parses a SMIRKS string as a reaction object and outputs a dictionary with 
    atom type and integer diff of product count - substrate count for each atom type,
    excluding zero differences.

    :param smirks: SMIRKS string
    :return: Dictionary with atom type and integer diff
    """
    rxn = AllChem.ReactionFromSmarts(smirks)
    atom_diff = {}

    # Count atoms in reactants
    reactant_atoms = {}
    for reactant in rxn.GetReactants():
        for atom in reactant.GetAtoms():
            atom_type = atom.GetSymbol()
            reactant_atoms[atom_type] = reactant_atoms.get(atom_type, 0) + 1

    # Count atoms in products
    product_atoms = {}
    for product in rxn.GetProducts():
        for atom in product.GetAtoms():
            atom_type = atom.GetSymbol()
            product_atoms[atom_type] = product_atoms.get(atom_type, 0) + 1

    # Calculate differences and exclude zeros
    all_atoms = set(reactant_atoms.keys()).union(set(product_atoms.keys()))
    for atom in all_atoms:
        diff = product_atoms.get(atom, 0) - reactant_atoms.get(atom, 0)
        if diff != 0:
            atom_diff[atom] = diff

    return atom_diff

def compare_atom_diffs(diff1: Dict[str, int], diff2: Dict[str, int]) -> bool:
    """
    Compares two atom type diff dictionaries to determine if they are the same,
    ignoring zero values.

    :param diff1: First atom type diff dictionary
    :param diff2: Second atom type diff dictionary
    :return: True if the dictionaries are the same, False otherwise
    """
    filtered_diff1 = {k: v for k, v in diff1.items() if v != 0}
    filtered_diff2 = {k: v for k, v in diff2.items() if v != 0}
    return filtered_diff1 == filtered_diff2

class AtomDiffCache:
    def __init__(self):
        self.csv_path = 'data/evodex_atom_diffs.csv'
        self.cache = {}
        self.load_csv()

    def load_csv(self):
        """
        Loads the CSV file and populates the cache with the dictionary
        where the key is the hash of the atom diff dictionary and the value is the operator name.
        """
        with open(self.csv_path, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                atom_diff = eval(row['atom_diff'])
                operator_name = row['operator_name']
                hash_key = self.hash_atom_diff(atom_diff)
                self.cache[hash_key] = operator_name

    @staticmethod
    def hash_atom_diff(atom_diff: Dict[str, int]) -> str:
        """
        Hashes the atom diff dictionary to create a unique key.

        :param atom_diff: Atom diff dictionary
        :return: Hash key as a string
        """
        atom_diff_str = str(sorted(atom_diff.items()))
        return hashlib.md5(atom_diff_str.encode()).hexdigest()

    def get_operator_name(self, atom_diff: Dict[str, int]) -> Optional[str]:
        """
        Retrieves the operator name for the given atom diff dictionary.

        :param atom_diff: Atom diff dictionary
        :return: Operator name if exists, None otherwise
        """
        hash_key = self.hash_atom_diff(atom_diff)
        return self.cache.get(hash_key)

# Tests for comparison logic
def test_compare_atom_diffs():
    tests = [
        # Test identical dictionaries
        ({'C': 1, 'O': -1}, {'C': 1, 'O': -1}, True),
        # Test differing values
        ({'C': 1, 'O': -1}, {'C': 1, 'O': -2}, False),
        # Test differing keys
        ({'C': 1, 'O': -1}, {'C': 1, 'N': -1}, False),
        # Test ignoring zero values
        ({'C': 1, 'O': -1}, {'C': 1, 'O': -1, 'H': 0}, True),
        ({'C': 0, 'O': 0}, {'C': 0, 'O': 0}, True),
        ({'C': 1}, {'C': 1, 'O': 0}, True),
        # Test empty dictionaries
        ({}, {}, True),
        # Test one empty, one non-empty
        ({}, {'C': 1}, False),
        # Test with zero values only
        ({'C': 0, 'O': 0}, {}, True),
    ]

    for i, (diff1, diff2, expected) in enumerate(tests):
        result = compare_atom_diffs(diff1, diff2)
        assert result == expected, f"Test {i + 1} failed: {diff1} vs {diff2} -> {result}, expected {expected}"
