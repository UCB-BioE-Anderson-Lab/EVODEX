import pytest
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os

# This deals with path issues
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from evodex.astatine import hydrogen_to_astatine_reaction, astatine_to_hydrogen_reaction

def get_carbon_count(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

@pytest.fixture(scope="module")
def load_reactions_with_hydrogen():
    reactions = []
    with open('tests/data/astatine_test_data.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append(row['rxn'])
    return reactions

def test_hydrogen_to_astatine_conversion(load_reactions_with_hydrogen):
    for reaction_smiles in load_reactions_with_hydrogen:
        mapped_smiles = reaction_smiles
            
        # Convert hydrogen to astatine
        astatine_smiles = hydrogen_to_astatine_reaction(mapped_smiles)
            
        # Convert astatine back to hydrogen
        reconstituted_smiles = astatine_to_hydrogen_reaction(astatine_smiles)
            
        # Check carbon counts in reactants and products
        reaction_orig = AllChem.ReactionFromSmarts(mapped_smiles, useSmiles=True)
        reaction_reconstituted = AllChem.ReactionFromSmarts(reconstituted_smiles, useSmiles=True)
            
        for orig_mol, recon_mol in zip(reaction_orig.GetReactants(), reaction_reconstituted.GetReactants()):
            orig_carbon_count = get_carbon_count(Chem.MolToSmiles(orig_mol, isomericSmiles=True))
            recon_carbon_count = get_carbon_count(Chem.MolToSmiles(recon_mol, isomericSmiles=True))
            assert orig_carbon_count == recon_carbon_count, f"Reactant carbon count mismatch: {orig_carbon_count} != {recon_carbon_count}"

        for orig_mol, recon_mol in zip(reaction_orig.GetProducts(), reaction_reconstituted.GetProducts()):
            orig_carbon_count = get_carbon_count(Chem.MolToSmiles(orig_mol, isomericSmiles=True))
            recon_carbon_count = get_carbon_count(Chem.MolToSmiles(recon_mol, isomericSmiles=True))
            assert orig_carbon_count == recon_carbon_count, f"Product carbon count mismatch: {orig_carbon_count} != {recon_carbon_count}"

if __name__ == "__main__":
    pytest.main()
