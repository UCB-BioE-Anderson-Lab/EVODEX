from rdkit import Chem
from rdkit.Chem import AllChem

def hydrogen_to_astatine_molecule(mol: Chem.Mol) -> Chem.Mol:
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:  # Hydrogen
            atom.SetAtomicNum(85)  # Astatine
    return mol

def astatine_to_hydrogen_molecule(mol: Chem.Mol) -> Chem.Mol:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 85:  # Astatine
            atom.SetAtomicNum(1)  # Hydrogen
    return mol

def hydrogen_to_astatine_reaction(reaction_smiles: str) -> str:
    reaction = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
    reactant_smiles = []
    product_smiles = []

    for mol in reaction.GetReactants():
        mol = hydrogen_to_astatine_molecule(mol)
        reactant_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    for mol in reaction.GetProducts():
        mol = hydrogen_to_astatine_molecule(mol)
        product_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    return '.'.join(reactant_smiles) + ">>" + '.'.join(product_smiles)

def astatine_to_hydrogen_reaction(reaction_smiles: str) -> str:
    reaction = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
    reactant_smiles = []
    product_smiles = []

    for mol in reaction.GetReactants():
        mol = astatine_to_hydrogen_molecule(mol)
        reactant_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    for mol in reaction.GetProducts():
        mol = astatine_to_hydrogen_molecule(mol)
        product_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    return '.'.join(reactant_smiles) + ">>" + '.'.join(product_smiles)

# # Test functions
# if __name__ == "__main__":
#     # Test on molecule
#     hydrogen_molecule = Chem.MolFromSmiles("CCO")
#     print("Hydrogen molecule SMILES:", Chem.MolToSmiles(hydrogen_molecule))

#     astatine_molecule = hydrogen_to_astatine_molecule(hydrogen_molecule)
#     print("Astatine molecule SMILES:", Chem.MolToSmiles(astatine_molecule))

#     reconstituted_hydrogen_molecule = astatine_to_hydrogen_molecule(astatine_molecule)
#     print("Reconstituted Hydrogen molecule SMILES:", Chem.MolToSmiles(reconstituted_hydrogen_molecule))

#     # Test on reaction
#     hydrogen_reaction_smiles = "[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH2:6]>>[CH3:1][C:2](=[O:3])[OH:6].[CH3:5][OH:4]"
#     print("Hydrogen reaction SMILES:", hydrogen_reaction_smiles)

#     astatine_reaction_smiles = hydrogen_to_astatine_reaction(hydrogen_reaction_smiles)
#     print("Astatine reaction SMILES:", astatine_reaction_smiles)

#     reconstituted_hydrogen_reaction_smiles = astatine_to_hydrogen_reaction(astatine_reaction_smiles)
#     print("Reconstituted Hydrogen reaction SMILES:", reconstituted_hydrogen_reaction_smiles)
