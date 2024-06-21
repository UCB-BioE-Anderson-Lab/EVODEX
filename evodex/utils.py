from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import InchiReadWriteError

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

def calculate_inchi(mol: Chem.Mol) -> str:
    """
    Calculates the InChI string for a given molecule.

    Parameters:
    mol (Chem.Mol): The RDKit molecule object.

    Returns:
    str: The InChI string of the molecule.

    Raises:
    InchiReadWriteError: If there's an error or warning that is not acceptable, with SMILES included.
    """
    try:
        inchi = Chem.MolToInchi(mol, treatWarningAsError=True)
        return inchi
    except InchiReadWriteError as e:
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception as smiles_error:
            smiles = f"Error generating SMILES: {smiles_error}"
        error_message = str(e)
        if any(warning in error_message for warning in ["Proton(s) added/removed", "Charges were rearranged"]):
            inchi = Chem.MolToInchi(mol, treatWarningAsError=False)
            return inchi
        raise InchiReadWriteError(f"InChI generation failed for molecule with SMILES: {smiles}. Error: {error_message}")