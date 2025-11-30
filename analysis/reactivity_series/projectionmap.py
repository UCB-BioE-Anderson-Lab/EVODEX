from rdkit import Chem
from rdkit.Chem import rdChemReactions
from evodex.astatine import hydrogen_to_astatine_molecule, astatine_to_hydrogen_molecule


def add_sequential_isotopes(smiles: str) -> str:
    """Add sequential isotopes starting at 1 to every atom in a SMILES string.

    Args:
        smiles: Input SMILES string for the substrate.

    Returns:
        A SMILES string with isotopes assigned sequentially.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")

    # Convert all hydrogens to astatine using the evodex.astatine utilities.
    # This also makes hydrogens explicit (via AddHs inside hydrogen_to_astatine_molecule),
    # so the subsequent isotope labeling can treat every atom uniformly.
    mol = hydrogen_to_astatine_molecule(mol, which_mol="substrate")

    for idx, atom in enumerate(mol.GetAtoms(), start=1):
        atom.SetIsotope(idx)

    # canonical SMILES with isotopes included
    isotopic_smiles = Chem.MolToSmiles(mol)
    return isotopic_smiles



def project_operator_to_mapped_products(substrate_smiles: str, operator_smirks: str) -> tuple[str, str]:
    """Project a reaction operator onto a substrate using isotopes as temporary labels,
    then transfer those isotopes to atom maps on both substrate and product.

    Returns:
        mapped_substrate_smiles, mapped_product_smiles
    """
    # First, label the substrate with sequential isotopes
    isotopic_substrate = add_sequential_isotopes(substrate_smiles)
    isotopic_sub_mol = Chem.MolFromSmiles(isotopic_substrate)
    if isotopic_sub_mol is None:
        raise ValueError(f"Could not parse isotopic substrate SMILES: {isotopic_substrate}")

    # Build the reaction from the SMIRKS
    rxn = rdChemReactions.ReactionFromSmarts(operator_smirks)
    if rxn is None:
        raise ValueError(f"Could not parse reaction operator SMIRKS: {operator_smirks}")

    # Apply the reaction to the isotopically labeled substrate
    products = rxn.RunReactants((isotopic_sub_mol,))
    if not products or not products[0]:
        raise ValueError("Reaction operator did not match isotopic substrate")

    product_mol = products[0][0]

    # Transfer isotopes to atom maps and clear isotopes on both substrate and product
    for mol in (isotopic_sub_mol, product_mol):
        for atom in mol.GetAtoms():
            iso = atom.GetIsotope()
            if iso:
                atom.SetAtomMapNum(iso)
                atom.SetIsotope(0)

    # Remove atom maps that do not appear on both sides
    sub_maps = {atom.GetAtomMapNum() for atom in isotopic_sub_mol.GetAtoms() if atom.GetAtomMapNum()}
    prod_maps = {atom.GetAtomMapNum() for atom in product_mol.GetAtoms() if atom.GetAtomMapNum()}
    non_matching = (sub_maps ^ prod_maps)

    if non_matching:
        for mol in (isotopic_sub_mol, product_mol):
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() in non_matching:
                    atom.SetAtomMapNum(0)

    # Convert astatine back to hydrogen on both substrate and product
    isotopic_sub_mol = astatine_to_hydrogen_molecule(isotopic_sub_mol)
    product_mol = astatine_to_hydrogen_molecule(product_mol)

    mapped_substrate_smiles = Chem.MolToSmiles(isotopic_sub_mol)
    mapped_product_smiles = Chem.MolToSmiles(product_mol)
    return mapped_substrate_smiles, mapped_product_smiles


def main():
    # Example substrate and operator
    substrate_smiles = "CCOC"
    operator_smirks = "[C:1][O:2][C]>>[C:1][O:2]"

    print("Substrate (input):", substrate_smiles)
    print("Operator (SMIRKS):", operator_smirks)

    mapped_substrate_smiles, mapped_product_smiles = project_operator_to_mapped_products(
        substrate_smiles, operator_smirks
    )

    print("Mapped substrate SMILES:", mapped_substrate_smiles)
    print("Mapped product SMILES:", mapped_product_smiles)


if __name__ == "__main__":
    main()