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

# def get_molecule_hash(mol: Chem.Mol) -> str:
#     """
#     Generates a hash for the given molecule after kekulization and InChI key generation.

#     Parameters:
#     mol (Chem.Mol): The RDKit molecule object.

#     Returns:
#     str: The InChI key of the molecule.
#     """
#     try:
#         mol = Chem.Mol(mol)  # Create a copy of the molecule
#         Chem.SanitizeMol(mol)  # Sanitize the molecule to calculate explicit valence and other properties
#         Chem.Kekulize(mol, clearAromaticFlags=True)  # Ensure kekulization
#         inchi_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
#         return inchi_key
#     except Exception as e:
#         raise RuntimeError(f"Failed to generate molecule hash: {e}")
  
from rdkit import Chem

def _is_sp3(atom):
    """
    Check if an atom is SP3 hybridized (all single bonds).
    """
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return False
    return True

def get_molecule_hash(smiles: str) -> str:
    # Parse the input SMILES to a mol
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string provided")

    # Instantiate a new RWMol
    rwmol = Chem.RWMol()

    # Iterate through the mol and create atoms in the RWMol
    atom_map = {}
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if 1 < atomic_num < 85:  # Exclude H (1) and At (85)
            new_atom = Chem.Atom(atomic_num)
            if atomic_num == 6:  # Only apply labeling to carbon atoms
                is_sp3 = _is_sp3(atom)
                if not is_sp3:
                    new_atom.SetIsotope(0)
                else:
                    new_atom.SetIsotope(1)
            new_idx = rwmol.AddAtom(new_atom)
            atom_map[atom.GetIdx()] = new_idx

    # Iterate through the bonds in the mol and add them to the RWMol as single bonds
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in atom_map and end_idx in atom_map:
            rwmol.AddBond(atom_map[begin_idx], atom_map[end_idx], Chem.BondType.SINGLE)

    # Generate canonical SMILES from the RWMol
    smiles_canonical = Chem.MolToSmiles(rwmol, canonical=True, allHsExplicit=False)
    return smiles_canonical
