from rdkit import Chem
from rdkit.Chem import rdChemReactions
from evodex.astatine import hydrogen_to_astatine_molecule, astatine_to_hydrogen_molecule

"""
Projection map inputs a reaction operator (SMIRKS) and a substrate (SMILES) and outputs a
fully mapped reaction smiles for the substrate-to-product reaction SMILES. For this, it first
converts the hydrogens to astatines such that all atoms are heavy and explicit. It then assigns
sequential atom map values to each substrate atom, including the astatine/hydrogens.
This information as stored in the isotope field of each atom. It then projects the operator 
on the substrate to generate the product SMILES. This will also project the isotope fields such
that these are now paired between substrate and product except where atoms are removed or added
during the reaction. The isotope fields are transferred to the atom maps, and unmatched atom maps
are removed. The astatines are converted back to hydrogens and the fully mapped partial reaction
is returned.
"""

def add_sequential_isotopes(smiles: str) -> str:
    """Assign sequential isotopes starting at 1 across all atoms in the molecule.

    The substrate's hydrogens are first converted to astatine to ensure that all
    atoms are heavy and explicit. Isotopes applied after conversion therefore
    provide a uniform index covering every atom.

    Args:
        smiles: Substrate SMILES.

    Returns:
        SMILES containing sequential isotopes.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")

    mol = hydrogen_to_astatine_molecule(mol, which_mol="substrate")

    for idx, atom in enumerate(mol.GetAtoms(), start=1):
        atom.SetIsotope(idx)

    return Chem.MolToSmiles(mol)


def project_operator_to_mapped_products(
    substrate_smiles: str, operator_smirks: str
) -> tuple[str, str]:
    """Apply an operator SMIRKS to a substrate and transfer isotope labels to atom maps.

    The substrate is first isotopically labeled. A reaction is then run on this
    labeled substrate. The first product of the first product set is retained.
    Isotopes are promoted to atom map numbers on both molecules, then removed
    when they appear on only one side.

    Args:
        substrate_smiles: Input substrate.
        operator_smirks: Reaction operator encoded as SMIRKS.

    Returns:
        (mapped_substrate_smiles, mapped_product_smiles)
    """
    isotopic_substrate = add_sequential_isotopes(substrate_smiles)
    isotopic_sub_mol = Chem.MolFromSmiles(isotopic_substrate)
    if isotopic_sub_mol is None:
        raise ValueError(f"Invalid isotopic substrate: {isotopic_substrate}")

    rxn = rdChemReactions.ReactionFromSmarts(operator_smirks)
    if rxn is None:
        raise ValueError(f"Could not parse operator SMIRKS: {operator_smirks}")

    products = rxn.RunReactants((isotopic_sub_mol,))
    if not products or not products[0]:
        raise ValueError("Operator did not match substrate")

    product_mol = products[0][0]

    for mol in (isotopic_sub_mol, product_mol):
        for atom in mol.GetAtoms():
            iso = atom.GetIsotope()
            if iso:
                atom.SetAtomMapNum(iso)
                atom.SetIsotope(0)

    sub_maps = {a.GetAtomMapNum() for a in isotopic_sub_mol.GetAtoms() if a.GetAtomMapNum()}
    prod_maps = {a.GetAtomMapNum() for a in product_mol.GetAtoms() if a.GetAtomMapNum()}
    non_shared = sub_maps ^ prod_maps

    if non_shared:
        for mol in (isotopic_sub_mol, product_mol):
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() in non_shared:
                    atom.SetAtomMapNum(0)

    isotopic_sub_mol = astatine_to_hydrogen_molecule(isotopic_sub_mol)
    product_mol = astatine_to_hydrogen_molecule(product_mol)

    return Chem.MolToSmiles(isotopic_sub_mol), Chem.MolToSmiles(product_mol)


if __name__ == "__main__":
    # Minimal example
    substrate = "CCOC"
    operator = "[C:1][O:2][C]>>[C:1][O:2]"
    sub_out, prod_out = project_operator_to_mapped_products(substrate, operator)
    print(sub_out)
    print(prod_out)
