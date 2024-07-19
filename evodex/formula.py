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

    # Switch 'At' to 'H' in the dictionary
    if 'At' in atom_diff:
        atom_diff['H'] = atom_diff.pop('At')

    return atom_diff

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
