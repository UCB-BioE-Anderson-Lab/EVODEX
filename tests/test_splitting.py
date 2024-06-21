import pytest
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
from evodex.splitting import split_reaction

@pytest.fixture(scope="module")
def load_reactions_with_astatine():
    reactions = []
    with open('tests/data/splitting_test_data.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append(row['atom_mapped'])
    return reactions

def test_split_reaction(load_reactions_with_astatine):
    for reaction_smiles in load_reactions_with_astatine:
        try:
            partial_reactions = split_reaction(reaction_smiles)
            
            # Ensure the list is not empty
            assert len(partial_reactions) > 0, "No partial reactions generated"
            
            for partial in partial_reactions:
                try:
                    rxn = AllChem.ReactionFromSmarts(partial, useSmiles=True)
                    
                    # Check that each partial reaction can be parsed as a reaction SMILES
                    assert rxn is not None, f"Invalid reaction SMILES: {partial}"
                    
                    # Ensure all molecules contain at least one mapped atom on both sides
                    substrate_atom_maps = set()
                    for mol in rxn.GetReactants():
                        for atom in mol.GetAtoms():
                            if atom.GetAtomMapNum() > 0:
                                substrate_atom_maps.add(atom.GetAtomMapNum())
                    
                    product_atom_maps = set()
                    for mol in rxn.GetProducts():
                        for atom in mol.GetAtoms():
                            if atom.GetAtomMapNum() > 0:
                                product_atom_maps.add(atom.GetAtomMapNum())
                    
                    assert len(substrate_atom_maps) > 0, f"No mapped atoms in reactants: {partial}"
                    assert len(product_atom_maps) > 0, f"No mapped atoms in products: {partial}"
                    
                    # Ensure all atom maps are unique on either side
                    assert len(substrate_atom_maps) == sum(1 for mol in rxn.GetReactants() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0), f"Duplicate atom maps in reactants: {partial}"
                    assert len(product_atom_maps) == sum(1 for mol in rxn.GetProducts() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0), f"Duplicate atom maps in products: {partial}"
                    
                    # Ensure the set of maps present in the substrates matches those in the products
                    assert substrate_atom_maps == product_atom_maps, f"Mismatch between reactant and product atom maps: {partial}"
                
                except Exception as e:
                    pytest.fail(f"Failed to validate partial reaction: {partial}, Error: {e}")
        
        except Exception as e:
            pytest.fail(f"split_reaction raised an error: {e}")

if __name__ == "__main__":
    pytest.main()
