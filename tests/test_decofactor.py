import pytest
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from evodex.decofactor import remove_cofactors

@pytest.fixture(scope="module")
def load_reactions_with_astatine():
    reactions = []
    with open('tests/data/splitting_test_data.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append(row['atom_mapped'])
    return reactions

def test_remove_cofactors(load_reactions_with_astatine):
    for reaction_smiles in load_reactions_with_astatine:
        try:
            partial_reaction = remove_cofactors(reaction_smiles)
            
            # Ensure the partial reaction is not None or empty
            assert partial_reaction is not None, "No partial reactions generated"
            
            # Check if the partial reaction is ">>"
            if partial_reaction == ">>":
                continue  # Accept ">>" as valid behavior
            
            try:
                rxn = rdChemReactions.ReactionFromSmarts(partial_reaction, useSmiles=True)
                
                # Check that each partial reaction can be parsed as a reaction SMILES
                assert rxn is not None, f"Invalid reaction SMILES: {partial_reaction}"
                
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
                
                assert len(substrate_atom_maps) > 0, f"No mapped atoms in reactants: {partial_reaction}"
                assert len(product_atom_maps) > 0, f"No mapped atoms in products: {partial_reaction}"
                
                # Ensure all atom maps are unique on either side
                assert len(substrate_atom_maps) == sum(1 for mol in rxn.GetReactants() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0), f"Duplicate atom maps in reactants: {partial_reaction}"
                assert len(product_atom_maps) == sum(1 for mol in rxn.GetProducts() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0), f"Duplicate atom maps in products: {partial_reaction}"
                
                # Ensure the set of maps present in the substrates matches those in the products
                assert substrate_atom_maps == product_atom_maps, f"Mismatch between reactant and product atom maps: {partial_reaction}"
            
            except Exception as e:
                pytest.fail(f"Failed to validate partial reaction: {partial_reaction}, Error: {e}")
        
        except Exception as e:
            pytest.fail(f"remove_cofactors raised an error: {e}")

if __name__ == "__main__":
    pytest.main()
