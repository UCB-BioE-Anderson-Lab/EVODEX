from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions

def hydrogen_to_astatine(reaction_smiles: str) -> str:
    reaction = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
    reactant_smiles = []
    product_smiles = []

    for mol in reaction.GetReactants():
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Hydrogen
                atom.SetAtomicNum(85)  # Astatine
        reactant_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    for mol in reaction.GetProducts():
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Hydrogen
                atom.SetAtomicNum(85)  # Astatine
        product_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    return '.'.join(reactant_smiles) + ">>" + '.'.join(product_smiles)

def astatine_to_hydrogen(reaction_smiles: str) -> str:
    reaction = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
    reactant_smiles = []
    product_smiles = []

    for mol in reaction.GetReactants():
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 85:  # Astatine
                atom.SetAtomicNum(1)  # Hydrogen
        reactant_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    for mol in reaction.GetProducts():
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 85:  # Astatine
                atom.SetAtomicNum(1)  # Hydrogen
        product_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    return '.'.join(reactant_smiles) + ">>" + '.'.join(product_smiles)

# # Test on SMILES
# hydrogen_smiles = "[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH2:6]>>[CH3:1][C:2](=[O:3])[OH:6].[CH3:5][OH:4]"
# print("Hydrogen SMILES:", hydrogen_smiles)

# astatine_smiles = hydrogen_to_astatine(hydrogen_smiles)
# print("Astatine SMILES:", astatine_smiles)

# reconstituted_hydrogen_smiles = astatine_to_hydrogen(astatine_smiles)
# print("Reconstituted Hydrogen SMILES:", reconstituted_hydrogen_smiles)

# # Test on SMIRKS
# hydrogen_smirks = "[H][C:1]([H:2])[O:3][H]>>[C:1]([H:2])=[O:3]"
# print("Hydrogen SMILES:", hydrogen_smirks)

# astatine_smirks = hydrogen_to_astatine(hydrogen_smirks)
# print("Astatine SMIRKS:", astatine_smirks)

# reconstituted_hydrogen_smiles = astatine_to_hydrogen(astatine_smirks)
# print("Reconstituted Hydrogen SMIRKS:", reconstituted_hydrogen_smiles)
