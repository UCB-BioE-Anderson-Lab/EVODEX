import csv
import hashlib
from typing import Dict, Optional
from rdkit import Chem
from rdkit.Chem import AllChem

# Exact mass values for elements
EXACT_MASS = {
    'H': 1.00783,
    'D': 2.01410,
    'C': 12.0000,
    'N': 14.0031,
    'O': 15.9949,
    'F': 18.9984,
    'Si': 27.9769,
    'P': 30.9738,
    'S': 31.9721,
    'Cl': 34.9689,
    'Br': 78.9183,
    'I': 126.9045
}

def calculate_formula_diff(smirks: str) -> Dict[str, int]:
    """
    Calculates the atom type difference for a SMIRKS string as a reaction object
    and outputs a dictionary with atom type and integer diff of product count - substrate count 
    for each atom type, excluding zero differences.

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

def _compare_atom_diffs(diff1: Dict[str, int], diff2: Dict[str, int]) -> bool:
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

def calculate_exact_mass(atom_diff: Dict[str, int]) -> float:
    """
    Calculates the exact mass for the given atom diff dictionary.

    :param atom_diff: Atom diff dictionary
    :return: Exact mass as a float
    """
    exact_mass = 0.0
    for atom, count in atom_diff.items():
        exact_mass += EXACT_MASS.get(atom, 0.0) * count
    return exact_mass

class AtomDiffCache:
    def __init__(self):
        self.csv_path = 'data/evodex_atom_diffs.csv'
        self.cache = {}
        self.mass_cache = {}
        self._load_csv()

    def _load_csv(self):
        """
        Loads the CSV file and populates the cache with the dictionary
        where the key is the hash of the atom diff dictionary and the value is the operator name.
        Also populates the mass_cache with exact masses for quick lookup.
        """
        with open(self.csv_path, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                atom_diff = eval(row['atom_diff'])
                operator_name = row['operator_name']
                hash_key = self._hash_atom_diff(atom_diff)
                self.cache[hash_key] = operator_name

                # Add exact mass to mass_cache
                exact_mass = calculate_exact_mass(atom_diff)
                self.mass_cache[exact_mass] = operator_name

    @staticmethod
    def _hash_atom_diff(atom_diff: Dict[str, int]) -> str:
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
        hash_key = self._hash_atom_diff(atom_diff)
        return self.cache.get(hash_key)

    def get_operator_by_mass(self, mass: float, resolution: float = 0.01) -> Optional[str]:
        """
        Retrieves the operator name for the given exact mass within the specified resolution.

        :param mass: Exact mass difference
        :param resolution: Resolution value in daltons
        :return: Operator name if exists within the resolution, None otherwise
        """
        for exact_mass, operator_name in self.mass_cache.items():
            if abs(exact_mass - mass) <= resolution:
                return operator_name
        return None
