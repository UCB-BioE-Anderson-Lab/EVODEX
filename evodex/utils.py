from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
import hashlib


def _is_sp3(atom):
    """
    Check if an atom is SP3 hybridized (all single bonds).
    """
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return False
    return True


def get_molecule_hash_from_mol(mol: Chem.Mol) -> str:
    """
    Generate a canonical, abstraction-stable hash string for an RDKit Mol.

    Behavior:
    - Drops H (1) and At (85)
    - Collapses all bonds to SINGLE
    - Encodes "oxidized" carbons (any non-single bond in input mol)
      via isotope: C-0 (oxidized) vs C-1 (sp3)
    - Does NOT parse or reparse SMILES
    """

    if mol is None:
        raise ValueError("Mol is None")

    rwmol = Chem.RWMol()
    atom_map = {}

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()

        # Exclude H (1), At (85), and any other oddities
        if atomic_num <= 1 or atomic_num >= 85:
            continue

        new_atom = Chem.Atom(atomic_num)

        oxidized = False
        for bond in atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                oxidized = True
                break

        if atomic_num == 6:
            new_atom.SetIsotope(0 if oxidized else 1)

        new_idx = rwmol.AddAtom(new_atom)
        atom_map[atom.GetIdx()] = new_idx

    for bond in mol.GetBonds():
        b = bond.GetBeginAtomIdx()
        e = bond.GetEndAtomIdx()
        if b in atom_map and e in atom_map:
            rwmol.AddBond(atom_map[b], atom_map[e], Chem.BondType.SINGLE)

    m2 = rwmol.GetMol()
    Chem.SanitizeMol(m2)

    return Chem.MolToSmiles(
        m2,
        canonical=True,
        allHsExplicit=False,
        isomericSmiles=False,
    )


def reaction_hash(smirks: str) -> int:
    """
    Compute a deterministic hash for a reaction SMIRKS.

    Key property:
    - Never treats SMARTS / SMIRKS as SMILES
    - Hashes directly from RDKit reaction templates
    """

    try:
        rxn = AllChem.ReactionFromSmarts(smirks)
        rxn.Initialize()
    except Exception as e:
        raise RuntimeError(
            f"Failed to load reaction from SMIRKS:\n{smirks}\n{e}"
        )

    def get_hash(mol: Chem.Mol) -> str:
        try:
            return get_molecule_hash_from_mol(mol)
        except Exception as e:
            raise RuntimeError(f"Failed to get molecule hash:\n{e}")

    try:
        substrate_hashes = {get_hash(m) for m in rxn.GetReactants()}
        product_hashes = {get_hash(m) for m in rxn.GetProducts()}
    except Exception as e:
        raise RuntimeError(
            f"Failed to generate hashes for substrates or products: {e}"
        )

    substrate_hash_str = "".join(sorted(substrate_hashes))
    product_hash_str = "".join(sorted(product_hashes))
    combined_hash_str = substrate_hash_str + ">>" + product_hash_str

    reaction_hash_value = hashlib.sha256(
        combined_hash_str.encode("utf-8")
    ).hexdigest()

    return int(reaction_hash_value, 16)