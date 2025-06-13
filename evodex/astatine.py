from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def copy_molecule_with_substitution(mol: Chem.Mol, substitution: dict, which_mol="") -> Chem.Mol:
    # Capture chiral tags BEFORE AddHs
    chiral_tags = {}
    for atom in mol.GetAtoms():
        atom_map_num = atom.GetAtomMapNum()
        if atom_map_num > 0:
            chiral_tags[atom_map_num] = atom.GetChiralTag()
    mol = Chem.AddHs(mol)

    new_mol = Chem.RWMol()
    atom_map = {}

    # Add atoms
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()

        if atomic_num == 1:
            # If substituting H → At
            if substitution.get(1, 1) == 85:
                new_atom = Chem.Atom(85)
                new_idx = new_mol.AddAtom(new_atom)
                atom_map[atom.GetIdx()] = new_idx

                new_mol.GetAtomWithIdx(new_idx).SetAtomMapNum(atom.GetAtomMapNum())
            # Else skip H atom (for At → H roundtrip)
            continue
        else:
            if atomic_num in substitution:
                atomic_num = substitution[atomic_num]

            new_atom = Chem.Atom(atomic_num)
            new_idx = new_mol.AddAtom(new_atom)
            atom_map[atom.GetIdx()] = new_idx

            new_mol.GetAtomWithIdx(new_idx).SetChiralTag(atom.GetChiralTag())
            new_mol.GetAtomWithIdx(new_idx).SetFormalCharge(atom.GetFormalCharge())
            new_mol.GetAtomWithIdx(new_idx).SetAtomMapNum(atom.GetAtomMapNum())

    # Add bonds
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        if begin_idx not in atom_map or end_idx not in atom_map:
            continue

        new_begin = atom_map[begin_idx]
        new_end = atom_map[end_idx]

        new_mol.AddBond(new_begin, new_end, bond.GetBondType())

    # Restore original chiral tags using atom map numbers
    for old_idx, new_idx in atom_map.items():
        atom_map_num = new_mol.GetAtomWithIdx(new_idx).GetAtomMapNum()
        if atom_map_num in chiral_tags:
            new_mol.GetAtomWithIdx(new_idx).SetChiralTag(chiral_tags[atom_map_num])

    result = new_mol.GetMol()
    Chem.SanitizeMol(result)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)

    return result

def hydrogen_to_astatine_molecule(mol: Chem.Mol, which_mol="") -> Chem.Mol:
    return copy_molecule_with_substitution(mol, {1: 85}, which_mol)

def astatine_to_hydrogen_molecule(mol: Chem.Mol) -> Chem.Mol:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 85:  # Astatine
            atom.SetAtomicNum(1)  # Hydrogen
    return mol

def hydrogen_to_astatine_reaction(reaction_smiles: str) -> str:
    reaction = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)

    reactant_smiles = []
    product_smiles = []

    for i, mol in enumerate(reaction.GetReactants()):
        mol = hydrogen_to_astatine_molecule(mol, which_mol=f"Reactant {i+1}")
        reactant_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True, canonical=False))

    for i, mol in enumerate(reaction.GetProducts()):
        mol = hydrogen_to_astatine_molecule(mol, which_mol=f"Product {i+1}")
        product_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True, canonical=False))

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

# DataFrame utility for converting a column of reaction SMILES using astatine_to_hydrogen_reaction
def convert_dataframe_smiles_column_at_to_h(df, column_name):
    """
    Applies astatine_to_hydrogen_reaction to the specified column of a DataFrame.
    Returns the converted DataFrame and a list of errors (tuples of index, value, and exception).
    """
    import pandas as pd
    errors = []
    def convert(row):
        try:
            return astatine_to_hydrogen_reaction(row)
        except Exception as e:
            errors.append((row, e))
            return None
    # Copy DataFrame to avoid modifying in place
    df_converted = df.copy()
    df_converted[column_name] = df_converted[column_name].apply(convert)
    return df_converted, errors

# DataFrame utility for converting a column of reaction SMIRKS using astatine_to_hydrogen_reaction
def convert_dataframe_smirks_column_at_to_h(df, column_name):
    """
    Applies astatine_to_hydrogen_molecule to each molecule in a SMIRKS reaction in the specified DataFrame column.
    Returns the converted DataFrame and a list of errors.
    """
    from rdkit.Chem import rdChemReactions
    errors = []

    def convert(row):
        try:
            reactionOperator = rdChemReactions.ReactionFromSmarts(row)
            new_reaction = rdChemReactions.ChemicalReaction()

            def convert_mol(mol):
                rw_mol = Chem.RWMol()
                atom_map = {}

                for atom in mol.GetAtoms():
                    atomic_num = atom.GetAtomicNum()
                    new_atom = Chem.Atom(1 if atomic_num == 85 else atomic_num)
                    new_idx = rw_mol.AddAtom(new_atom)
                    atom_map[atom.GetIdx()] = new_idx

                    rw_atom = rw_mol.GetAtomWithIdx(new_idx)
                    rw_atom.SetChiralTag(atom.GetChiralTag())
                    rw_atom.SetFormalCharge(atom.GetFormalCharge())
                    rw_atom.SetAtomMapNum(atom.GetAtomMapNum())

                for bond in mol.GetBonds():
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()

                    if begin_idx not in atom_map or end_idx not in atom_map:
                        continue

                    new_begin = atom_map[begin_idx]
                    new_end = atom_map[end_idx]

                    rw_mol.AddBond(new_begin, new_end, bond.GetBondType())

                result = rw_mol.GetMol()
                Chem.SanitizeMol(result)
                return result

            for i in range(reactionOperator.GetNumReactantTemplates()):
                mol = reactionOperator.GetReactantTemplate(i)
                converted = convert_mol(mol)
                new_reaction.AddReactantTemplate(converted)

            for i in range(reactionOperator.GetNumProductTemplates()):
                mol = reactionOperator.GetProductTemplate(i)
                converted = convert_mol(mol)
                new_reaction.AddProductTemplate(converted)

            smirks = rdChemReactions.ReactionToSmarts(new_reaction)
            return smirks
        except Exception as e:
            errors.append((row, e))
            return None

    df_converted = df.copy()
    df_converted[column_name] = df_converted[column_name].apply(convert)
    return df_converted, errors


# Allow running and testing directly
if __name__ == "__main__":
    # Example test reaction (replace with your test case)
    test_reaction = "[CH3:1][O:2][C:3](=[O:4])[C@H:5]([CH2:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1)[NH:13][C:14](=[O:15])[C@H:16]([CH2:17][C:18](=[O:19])[OH:20])[NH:21][C:22](=[O:23])[O:24][CH2:25][c:26]1[cH:27][cH:28][cH:29][cH:30][cH:31]1.[OH2:32]>>[C:14](=[O:15])([C@H:16]([CH2:17][C:18](=[O:19])[OH:20])[NH:21][C:22](=[O:23])[O:24][CH2:25][c:26]1[cH:27][cH:28][cH:29][cH:30][cH:31]1)[OH:32].[CH3:1][O:2][C:3](=[O:4])[C@H:5]([CH2:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1)[NH2:13]"

    print("Original reactionpython -m pipeline.phase6_synthesis_subset:")
    print(test_reaction)

    astatine_rxn = hydrogen_to_astatine_reaction(test_reaction)
    print("\nHydrogen to Astatine:")
    print(astatine_rxn)

    roundtrip_rxn = astatine_to_hydrogen_reaction(astatine_rxn)
    print("\nRoundtrip back to Hydrogen:")
    print(roundtrip_rxn)
