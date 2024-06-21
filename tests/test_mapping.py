import pytest
import csv
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from evodex.mapping import map_atoms

@pytest.fixture(scope="module")
def load_reactions_with_astatine():
    reactions = []
    with open('tests/data/mapping_test_data.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append(row['astatine_mapped'])
    return reactions

def get_atom_maps(mol):
    atom_maps = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 85:  # Exclude astatine
            atom_map = atom.GetAtomMapNum()
            if atom_map > 0:
                atom_maps.add(atom_map)
    return atom_maps

def test_map_atoms(load_reactions_with_astatine):
    for reaction_smiles in load_reactions_with_astatine:
        try:
            result = map_atoms(reaction_smiles)
            # Ensure Astatine is present in the result
            assert "At" in result  
            # Ensure it's a reaction SMILES
            assert ">>" in result  

            # Parse the reaction and validate atom maps
            reaction = rdChemReactions.ReactionFromSmarts(result, useSmiles=True)
            reactant_maps = set()
            product_maps = set()

            for reactant in reaction.GetReactants():
                reactant_maps.update(get_atom_maps(reactant))
            
            for product in reaction.GetProducts():
                product_maps.update(get_atom_maps(product))

            # Ensure all non-At atoms are mapped
            assert len(reactant_maps) > 0, "No reactant atom maps found"
            assert len(product_maps) > 0, "No product atom maps found"

            # Ensure atom maps are consistent between reactants and products
            assert reactant_maps == product_maps, f"Mismatch between reactant and product atom maps in {reaction_smiles}"

            # Ensure atom maps are unique within reactants and products
            for reactant in reaction.GetReactants():
                atom_maps = get_atom_maps(reactant)
                assert len(atom_maps) == len(set(atom_maps)), f"Duplicate atom maps in reactant {reaction_smiles}"
            
            for product in reaction.GetProducts():
                atom_maps = get_atom_maps(product)
                assert len(atom_maps) == len(set(atom_maps)), f"Duplicate atom maps in product {reaction_smiles}"

        except ValueError as e:
            pytest.fail(f"map_atoms raised an error: {e}")

if __name__ == "__main__":
    pytest.main()
