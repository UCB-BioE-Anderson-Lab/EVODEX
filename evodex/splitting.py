from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from itertools import combinations

def reaction_hash(rxn: AllChem.ChemicalReaction) -> tuple:
    substrate_inchis = set(Chem.MolToInchi(mol) for mol in rxn.GetReactants())
    product_inchis = set(Chem.MolToInchi(mol) for mol in rxn.GetProducts())
    return (frozenset(substrate_inchis), frozenset(product_inchis))

def split_reaction(smirks: str) -> list[str]:
    # Load the input SMIRKS as a reaction object
    rxn = AllChem.ReactionFromSmarts(smirks)
    
    # Standardize aromaticity and calculate reaction centers
    for mol in rxn.GetReactants():
        Chem.SanitizeMol(mol)
    for mol in rxn.GetProducts():
        Chem.SanitizeMol(mol)
    
    # Get the reactants and products as lists
    reactants = [Chem.MolToSmiles(mol) for mol in rxn.GetReactants()]
    products = [Chem.MolToSmiles(mol) for mol in rxn.GetProducts()]
    
    # Construct new reaction objects combinatorially
    reaction_splits = []
    for i in range(1, len(reactants) + 1):
        for j in range(1, len(products) + 1):
            for reactant_combo in combinations(reactants, i):
                for product_combo in combinations(products, j):
                    reactant_smiles = '.'.join(reactant_combo)
                    product_smiles = '.'.join(product_combo)
                    reaction_splits.append(f"{reactant_smiles}>>{product_smiles}")

    # Prune out reactions with no matching atom maps
    pruned_reactions = []
    for reaction in reaction_splits:
        try:
            rxn = AllChem.ReactionFromSmarts(reaction)
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
            
            if good_reaction:
                pruned_reactions.append(rxn)
        except Exception as e:
            print(f"Failed to process reaction: {reaction}, Error: {e}")

    # Process pruned reactions to clean up atom maps
    cleaned_reactions = []
    for rxn in pruned_reactions:
        try:
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
            
            # Remove unmapped molecules
            reactants = [mol for mol in rxn.GetReactants() if any(atom.GetAtomMapNum() > 0 for atom in mol.GetAtoms())]
            products = [mol for mol in rxn.GetProducts() if any(atom.GetAtomMapNum() > 0 for atom in mol.GetAtoms())]
            
            if reactants and products:
                cleaned_rxn = AllChem.ChemicalReaction()
                for mol in reactants:
                    cleaned_rxn.AddReactantTemplate(mol)
                for mol in products:
                    cleaned_rxn.AddProductTemplate(mol)
                cleaned_reactions.append(cleaned_rxn)

        except Exception as e:
            print(f"Failed to clean reaction: {AllChem.ReactionToSmarts(rxn)}, Error: {e}")

    # Remove duplicate reactions
    unique_reactions = set()
    final_reactions = []
    for rxn in cleaned_reactions:
        rxn_hash = reaction_hash(rxn)
        if rxn_hash not in unique_reactions:
            unique_reactions.add(rxn_hash)
            final_reactions.append(AllChem.ReactionToSmarts(rxn))

    return final_reactions