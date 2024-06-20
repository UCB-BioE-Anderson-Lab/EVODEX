from rdkit import Chem
from rdkit.Chem import AllChem

def validate_smiles(smiles: str) -> bool:
    """
    Validates if the given SMILES string is a valid chemical reaction.

    Parameters:
    smiles (str): The SMILES string to validate.

    Returns:
    bool: True if the SMILES string is valid, False otherwise.
    """
    try:
        reaction = AllChem.ReactionFromSmarts(smiles, useSmiles=True)
        if reaction is None:
            return False
        return True
    except:
        return False
