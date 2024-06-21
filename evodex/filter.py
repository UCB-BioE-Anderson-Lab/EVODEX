from rdkit import Chem
from rdkit.Chem import AllChem
from evodex.utils import calculate_inchi

def validate_smiles(smiles: str) -> bool:
    """
    Check if the given SMILES string passes the InChI filtering test.

    Parameters:
    smiles (str): The SMILES string to validate.

    Returns:
    bool: True if the SMILES string passes the filtering test, False otherwise.
    """
    try:
        reaction = AllChem.ReactionFromSmarts(smiles, useSmiles=True)
        reactants = list(reaction.GetReactants())
        products = list(reaction.GetProducts())

        for mol in reactants + products:
            calculate_inchi(mol)
    except Exception as e:
        return False

    return True