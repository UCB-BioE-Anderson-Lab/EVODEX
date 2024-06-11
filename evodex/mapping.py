from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import copy

def map_hydrogens(smiles: str) -> str:
    """
    This function takes a SMILES string representing a chemical reaction,
    validates and modifies it by adding explicit hydrogens, and ensures that
    all atoms, including most hydrogens, have unique and valid atom maps.

    Steps:
    1. Parse the SMILES string into a reaction object.
    2. Add explicit hydrogens to all reactants and products.
    3. Create a new reaction object with the added hydrogens.
    4. Validate that all heavy atoms have unique and valid atom maps.
    5. Ensure that atom maps are consistent between reactants and products.
    6. Assign atom maps to hydrogens, maintaining consistency.
    7. Convert the modified reaction back to a SMILES string.
    8. Return the modified SMILES string along with the original and modified reaction objects.
    """

    # Read the reaction SMILES and convert to a mol object
    try:
        raw_rxn = AllChem.ReactionFromSmarts(smiles, useSmiles=True)  # as smiles
        # print("Reaction successfully parsed.")
    except Exception as e:
        raise ValueError(f"Invalid reaction SMILES: {smiles}") from e

    # Add explicit hydrogens to all reactants and products
    new_reactants = []
    for reactant in raw_rxn.GetReactants():
        reactant_with_h = Chem.AddHs(reactant)
        new_reactants.append(reactant_with_h)

    new_products = []
    for product in raw_rxn.GetProducts():
        product_with_h = Chem.AddHs(product)
        new_products.append(product_with_h)

    # Create a new reaction with added hydrogens
    reaction = AllChem.ChemicalReaction()
    for reactant in new_reactants:
        reaction.AddReactantTemplate(reactant)
    for product in new_products:
        reaction.AddProductTemplate(product)

    # Make a duplicate of the original for output
    original_reaction = copy.deepcopy(reaction)

    # Validate atom maps
    def validate_atom_maps(mol):
        atom_map_set = set()
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Skip hydrogen
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
        reactants_atom_maps.update(validate_atom_maps(reactant))
        # print(f"Reactants atom maps: {reactants_atom_maps}")

    for product in reaction.GetProducts():
        products_atom_maps.update(validate_atom_maps(product))
        # print(f"Products atom maps: {products_atom_maps}")

    if reactants_atom_maps != products_atom_maps:
        raise ValueError(f"Mismatch between reactant and product atom maps in {smiles}")

    # Add atom maps to hydrogens
    next_atom_map = max(reactants_atom_maps.union(products_atom_maps)) + 1
    hydrogen_map_dict = {}
    # print(f"Initial next available atom map: {next_atom_map}")

    def add_hydrogen_maps(mol, is_reactant=True):
        nonlocal next_atom_map
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Hydrogen
                for neighbor in atom.GetNeighbors():
                    neighbor_map = neighbor.GetAtomMapNum()
                    if neighbor_map > 0:
                        if is_reactant:
                            if neighbor_map not in hydrogen_map_dict:
                                hydrogen_map_dict[neighbor_map] = []
                            hydrogen_map_dict[neighbor_map].append(next_atom_map)
                            atom.SetAtomMapNum(next_atom_map)
                            # print(f"Assigned atom map {next_atom_map} to reactant hydrogen attached to atom map {neighbor_map}")
                            next_atom_map += 1
                        else:
                            if neighbor_map in hydrogen_map_dict and hydrogen_map_dict[neighbor_map]:
                                assigned_map = hydrogen_map_dict[neighbor_map].pop(0)
                                atom.SetAtomMapNum(assigned_map)
                                # print(f"Assigned atom map {assigned_map} to product hydrogen attached to atom map {neighbor_map}")
                            else:
                                atom.SetAtomMapNum(0)
                                # print(f"Assigned atom map 0 to product hydrogen attached to atom map {neighbor_map} (no matching reactant hydrogen)")

    for reactant in reaction.GetReactants():
        add_hydrogen_maps(reactant, is_reactant=True)

    for product in reaction.GetProducts():
        add_hydrogen_maps(product, is_reactant=False)

    for reactant in reaction.GetReactants():
        for atom in reactant.GetAtoms():
            if atom.GetAtomicNum() == 1 and atom.GetAtomMapNum() in [item for sublist in hydrogen_map_dict.values() for item in sublist]:
                atom.SetAtomMapNum(0)
                # print(f"Reverted reactant hydrogen atom map to 0 for atom index {atom.GetIdx()}")

    # Make a duplicate of the original for output
    h_mapped_reaction = copy.deepcopy(reaction)

    # Convert mol object back to SMIRKS string
    try:
        modified_smiles = AllChem.ReactionToSmarts(reaction)
    except Exception as e:
        raise ValueError(f"Error converting modified reaction to SMIRKS: {smiles}") from e

    return modified_smiles, original_reaction, h_mapped_reaction
