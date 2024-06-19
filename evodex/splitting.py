from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import combinations

def canonicalize_molecule(mol: Chem.Mol) -> str:
    try:
        mol_copy = Chem.Mol(mol)  # Make a copy of the molecule
        for atom in mol_copy.GetAtoms():  # Clear atom maps
            atom.SetAtomMapNum(0)
        for atom in mol_copy.GetAtoms():  # Update implicit valence
            atom.UpdatePropertyCache()
        mol_copy = Chem.AddHs(mol_copy)  # Add explicit hydrogens
        Chem.SanitizeMol(mol_copy)  # Sanitize the molecule
        canonical_smiles = Chem.MolToSmiles(mol_copy, canonical=True)
    except Exception as e:
        raise RuntimeError(f"Failed to canonicalize molecule: {e}")
    return canonical_smiles

def reaction_hash(rxn: AllChem.ChemicalReaction) -> tuple:
    try:
        substrate_smiles = set(canonicalize_molecule(mol) for mol in rxn.GetReactants())
        product_smiles = set(canonicalize_molecule(mol) for mol in rxn.GetProducts())
    except Exception as e:
        raise RuntimeError(f"Failed to generate reaction hash: {e}")
    return (frozenset(substrate_smiles), frozenset(product_smiles))

def split_reaction(smirks: str) -> list[str]:
    try:
        # Load the input SMIRKS as a reaction object
        rxn = AllChem.ReactionFromSmarts(smirks)
        if not rxn:
            raise ValueError(f"Invalid SMIRKS string: {smirks}")

        # Get the reactants and products as lists
        reactants = [Chem.MolToSmiles(mol, isomericSmiles=True) for mol in rxn.GetReactants()]
        products = [Chem.MolToSmiles(mol, isomericSmiles=True) for mol in rxn.GetProducts()]

        print("reactants:", reactants)
        print("products:", products)

        # Construct new reaction objects combinatorially
        reaction_splits = []
        for i in range(1, len(reactants) + 1):
            for j in range(1, len(products) + 1):
                for reactant_combo in combinations(reactants, i):
                    for product_combo in combinations(products, j):
                        reactant_smiles = '.'.join(reactant_combo)
                        product_smiles = '.'.join(product_combo)
                        reaction_splits.append(f"{reactant_smiles}>>{product_smiles}")

        print("reaction_splits:", len(reaction_splits))

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
                raise RuntimeError(f"Failed to process reaction: {reaction}, Error: {e}")

        print("pruned_reactions:", len(pruned_reactions))

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
                raise RuntimeError(f"Failed to clean reaction: {AllChem.ReactionToSmarts(rxn)}, Error: {e}")

        # Remove duplicate reactions
        unique_reactions = set()
        final_reactions = []
        for rxn in cleaned_reactions:
            rxn_hash = reaction_hash(rxn)
            if rxn_hash not in unique_reactions:
                unique_reactions.add(rxn_hash)
                final_reactions.append(AllChem.ReactionToSmarts(rxn))

    except Exception as e:
        raise RuntimeError(f"Failed to split reaction: {e}")

    return final_reactions
