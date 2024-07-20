from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import combinations
import logging

logging.basicConfig(level=logging.INFO)

def split_reaction(smiles: str) -> list[str]:
    """
    This script reads in a full reaction and generates EVODEX-P (partial reactions) from it.
    It does this by enumerating all combinations of substrates and products, and then filters
    out any partials that don't satisfy the criteria of having at least one matching atom map
    between reactants and products. The algorithm proceeds as follows:

    1. **Splitting SMILES**: The reaction SMILES string is split into substrates (reactants) and
       products. Each part is further split into individual molecules.
    2. **Index Assignment**: Integer indices are assigned to each reactant and product molecule.
    3. **Combination Construction**: All possible combinations of reactants and products are generated.
    4. **Pruning Reactions**: Reactions without matching atom maps between reactants and products
       are filtered out. This ensures only chemically meaningful reactions are retained.
    5. **Cleaning Atom Maps**: Atom maps are adjusted to remove unnecessary mappings, ensuring
       consistency in the representation.
    6. **Error Handling**: Logging is used to capture and report errors encountered during processing.

    The result is a list of SMILES strings representing the valid partial reactions.

    Parameters:
    smiles (str): A reaction SMILES string in the form "reactants>>products".

    Returns:
    list[str]: A list of SMILES strings representing the valid partial reactions.
    """

    # Split the SMILES string into substrates (reactants) and products
    substrates, products = smiles.split('>>')
    substrate_list = substrates.split('.')
    product_list = products.split('.')

    # Assign an integer index to each substrate and product
    substrate_indices = list(range(len(substrate_list)))
    product_indices = list(range(len(product_list)))

    # Construct new reaction objects combinatorially
    all_partials = set()
    for i in range(1, len(substrate_indices) + 1):
        for j in range(1, len(product_indices) + 1):
            for reactant_combo in combinations(substrate_indices, i):
                for product_combo in combinations(product_indices, j):
                    all_partials.add((frozenset(reactant_combo), frozenset(product_combo)))

    # Prune out reactions with no matching atom maps
    pruned_reactions = set()
    for partial in all_partials:
        reactant_indices, product_indices = partial
        reactant_smiles = '.'.join([substrate_list[i] for i in sorted(reactant_indices)])
        product_smiles = '.'.join([product_list[i] for i in sorted(product_indices)])
        reaction_smiles = f"{reactant_smiles}>>{product_smiles}"
        
        try:
            # Create a reaction object from the SMILES
            rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
            substrate_atom_maps = set()
            # Collect atom maps from reactants
            for mol in rxn.GetReactants():
                for atom in mol.GetAtoms():
                    atom_map_num = atom.GetAtomMapNum()
                    if atom_map_num > 0:
                        substrate_atom_maps.add(atom_map_num)

            # Check for matching atom maps in products
            good_reaction = False
            for mol in rxn.GetProducts():
                for atom in mol.GetAtoms():
                    if atom.GetAtomMapNum() in substrate_atom_maps:
                        good_reaction = True
                        break
                if good_reaction:
                    break

            # If a valid reaction is found, add it to pruned reactions
            if good_reaction:
                pruned_reactions.add((frozenset(reactant_indices), frozenset(product_indices)))
        except Exception as e:
            logging.error(f"Failed to process reaction: {reaction_smiles}, Error: {e}")

    # Process pruned reactions to clean up atom maps
    cleaned_smiles = []
    for partial in pruned_reactions:
        reactant_indices, product_indices = partial
        reactant_smiles = '.'.join([substrate_list[i] for i in sorted(reactant_indices)])
        product_smiles = '.'.join([product_list[i] for i in sorted(product_indices)])
        reaction_smiles = f"{reactant_smiles}>>{product_smiles}"
        
        try:
            # Create a reaction object from the SMILES
            rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
            substrate_atom_maps = set()

            # Collect atom maps from reactants
            for mol in rxn.GetReactants():
                for atom in mol.GetAtoms():
                    atom_map_num = atom.GetAtomMapNum()
                    if atom_map_num > 0:
                        substrate_atom_maps.add(atom_map_num)

            # Adjust atom maps in products
            for mol in rxn.GetProducts():
                for atom in mol.GetAtoms():
                    atom_map_num = atom.GetAtomMapNum()
                    if atom_map_num > 0:
                        if atom_map_num not in substrate_atom_maps:
                            atom.SetAtomMapNum(0)
                        else:
                            substrate_atom_maps.remove(atom_map_num)

            # Adjust atom maps in reactants
            for mol in rxn.GetReactants():
                for atom in mol.GetAtoms():
                    atom_map_num = atom.GetAtomMapNum()
                    if atom_map_num in substrate_atom_maps:
                        atom.SetAtomMapNum(0)

            # Check for unmapped molecules
            has_unmapped_molecules = (
                any(not any(atom.GetAtomMapNum() > 0 for atom in mol.GetAtoms()) for mol in rxn.GetReactants()) or
                any(not any(atom.GetAtomMapNum() > 0 for atom in mol.GetAtoms()) for mol in rxn.GetProducts())
            )

            if not has_unmapped_molecules:
                cleaned_rxn = AllChem.ChemicalReaction()
                for mol in rxn.GetReactants():
                    cleaned_rxn.AddReactantTemplate(mol)
                for mol in rxn.GetProducts():
                    cleaned_rxn.AddProductTemplate(mol)
                cleaned_smiles.append(AllChem.ReactionToSmarts(cleaned_rxn))
        except Exception as e:
            logging.error(f"Failed to clean reaction: {reaction_smiles}, Error: {e}")

    return cleaned_smiles
