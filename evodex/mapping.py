from rdkit import Chem
from rdkit.Chem import AllChem
import copy

def map_atoms(smiles: str) -> str:
    """
    This function takes a SMILES string representing a chemical reaction,
    modifies it by ensuring that all atoms, including astatine, have unique
    and valid atom maps.

    Parameters:
    smiles (str): The input SMILES string with astatines replacing hydrogen

    Returns:
    str: The modified SMILES string with atom maps.
    """
    try:
        raw_rxn = AllChem.ReactionFromSmarts(smiles, useSmiles=True)
    except Exception as e:
        raise ValueError(f"Invalid reaction SMILES: {smiles}") from e

    new_reactants = list(raw_rxn.GetReactants())
    new_products = list(raw_rxn.GetProducts())

    # Calculate InChI key for all reactants and products
    reactant_inchis = [Chem.MolToInchiKey(mol) for mol in new_reactants]
    product_inchis = [Chem.MolToInchiKey(mol) for mol in new_products]

    # Remove identical InChI pairs from reactants and products
    to_remove_reactants = []
    to_remove_products = []
    for i, reactant_inchi in enumerate(reactant_inchis):
        for j, product_inchi in enumerate(product_inchis):
            if reactant_inchi == product_inchi:
                to_remove_reactants.append(i)
                to_remove_products.append(j)
                break

    new_reactants = [mol for i, mol in enumerate(new_reactants) if i not in to_remove_reactants]
    new_products = [mol for j, mol in enumerate(new_products) if j not in to_remove_products]

    reaction = AllChem.ChemicalReaction()
    for reactant in new_reactants:
        reaction.AddReactantTemplate(reactant)
    for product in new_products:
        reaction.AddProductTemplate(product)

    original_reaction = copy.deepcopy(reaction)

    def _validate_atom_maps(mol):
        atom_map_set = set()
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 85:  # Astatine
                continue
            atom_map = atom.GetAtomMapNum()
            if atom_map <= 0:
                raise ValueError(f"Atom without valid map number: {atom.GetSymbol()} in {smiles}")
            if atom_map in atom_map_set:
                raise ValueError(f"Duplicate atom map number: {atom_map} in {smiles}")
            atom_map_set.add(atom_map)
        return atom_map_set

    reactants_atom_maps = set()
    products_atom_maps = set()

    for reactant in reaction.GetReactants():
        reactants_atom_maps.update(_validate_atom_maps(reactant))

    for product in reaction.GetProducts():
        products_atom_maps.update(_validate_atom_maps(product))

    if reactants_atom_maps != products_atom_maps:
        raise ValueError(f"Mismatch between reactant and product atom maps in {smiles}")

    next_atom_map = max(reactants_atom_maps.union(products_atom_maps)) + 1
    astatine_map_dict = {}

    def _add_astatine_maps(mol, is_reactant=True):
        nonlocal next_atom_map
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 85:  # Astatine
                for neighbor in atom.GetNeighbors():
                    neighbor_map = neighbor.GetAtomMapNum()
                    if neighbor_map > 0:
                        if is_reactant:
                            if neighbor_map not in astatine_map_dict:
                                astatine_map_dict[neighbor_map] = []
                            astatine_map_dict[neighbor_map].append(next_atom_map)
                            atom.SetAtomMapNum(next_atom_map)
                            next_atom_map += 1
                        else:
                            if neighbor_map in astatine_map_dict and astatine_map_dict[neighbor_map]:
                                assigned_map = astatine_map_dict[neighbor_map].pop(0)
                                atom.SetAtomMapNum(assigned_map)
                            else:
                                atom.SetAtomMapNum(0)

    for reactant in reaction.GetReactants():
        _add_astatine_maps(reactant, is_reactant=True)

    for product in reaction.GetProducts():
        _add_astatine_maps(product, is_reactant=False)

    for reactant in reaction.GetReactants():
        for atom in reactant.GetAtoms():
            if atom.GetAtomicNum() == 85 and atom.GetAtomMapNum() in [item for sublist in astatine_map_dict.values() for item in sublist]:
                atom.SetAtomMapNum(0)

    try:
        reaction_smarts = AllChem.ReactionToSmarts(reaction)
        reactant_smiles = [Chem.MolToSmarts(mol, isomericSmiles=True) for mol in reaction.GetReactants()]
        product_smiles = [Chem.MolToSmarts(mol, isomericSmiles=True) for mol in reaction.GetProducts()]
        modified_smiles = '>>'.join(['.'.join(reactant_smiles), '.'.join(product_smiles)])
    except Exception as e:
        raise ValueError(f"Error converting modified reaction to SMIRKS: {reaction_smarts}") from e

    return modified_smiles
