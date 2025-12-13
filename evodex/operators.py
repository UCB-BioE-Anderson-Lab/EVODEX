from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import RDLogger
import re

RDLogger.DisableLog("rdApp.*")


# ================================================================
# operator_extractor: extract an operator SMIRKS from a reaction SMIRKS
# ================================================================

def operator_extractor(
    smirks: str,
    include_stereochemistry: bool = True,
    include_sigma: bool = True,
    include_pi: bool = True,
    include_extended: bool = True,
    include_unmapped: bool = True,
):
    """
    Extract a mechanism-focused operator from a full reaction SMIRKS.

    Steps
    1. Prepare reaction: parse SMIRKS and initialize (stereochemistry handled at the end).
    2. Reactive centers (bonding): mapped atoms whose bonding changes between sides; helper also adds directly attached unmapped neighbors for boundary stability.
    3. Covalent shell (sigma): atoms one σ bond away from reactive centers.
    4. Delocalized shell (pi): growth across conjugated/unsaturated connectivity; mirrored via atom-map numbers; excludes inductive effects and sp3-only participation.
    5. Extended shell: σ neighbors of the delocalized shell, excluding atoms already in covalent/delocalized shells.
    6. Identify unmapped atoms: collect atoms lacking map numbers on both sides.
    7. Strip stereochemistry (optional): if disabled, remove all stereo from a copy of the reaction.
    8. Compile operator: assemble keep sets, prune others, and emit SMIRKS.
    """

    # Step 1: Prepare reaction (parse and initialize; do not strip stereochemistry here)
    reaction = rdChemReactions.ReactionFromSmarts(smirks)
    reaction.Initialize()

    # Step 2: Bond-Changing reactive centers (bonding; helper also adds σ-adjacent unmapped neighbors)
    reactive_centers = _compute_reactive_centers(reaction)

    # Step 3: Covalent shell (σ/covalent shell: radius-1 from reactive centers)
    if include_sigma:
        covalent_shell_indices = _compute_covalent_shells(reaction, reactive_centers)
    else:
        covalent_shell_indices = None

    # Step 4: Delocalized shell (π/delocalized shell: conjugation-driven; mirrored across map numbers)
    if include_sigma and include_pi:
        delocalized_shell_indices = _grow_delocalized_shell(
            reaction,
            covalent_shell_indices,
        )
    else:
        delocalized_shell_indices = None

    # Step 5: Extended shell (σ neighbors of delocalized shell, excluding existing covalent/delocalized members)
    if include_sigma and include_pi and include_extended:
        extended_shell_indices = _add_extended_shell(
            reaction,
            covalent_shell_indices,
            delocalized_shell_indices,
        )
    else:
        extended_shell_indices = None

    # Step 6: Identify unmapped atoms on reactant and product sides (atoms lacking map numbers)
    unmapped_indices = _compute_unmapped_sets(reaction)
    if not include_unmapped:
        unmapped_indices = None

    # Step 7: Strip stereochemistry if requested
    reaction_for_compile = reaction
    if not include_stereochemistry:
        reaction_for_compile = rdChemReactions.ChemicalReaction()
        for i in range(reaction.GetNumReactantTemplates()):
            r = Chem.Mol(reaction.GetReactantTemplate(i))
            r = _strip_stereochemistry(r)
            reaction_for_compile.AddReactantTemplate(r)
        for i in range(reaction.GetNumProductTemplates()):
            p = Chem.Mol(reaction.GetProductTemplate(i))
            p = _strip_stereochemistry(p)
            reaction_for_compile.AddProductTemplate(p)
        reaction_for_compile.Initialize()

    # Step 8: Compile operator
    out_smirks = _compile_operator(
        reaction=reaction_for_compile,
        reactive_centers=reactive_centers,
        covalent_shell_indices=covalent_shell_indices,
        delocalized_shell_indices=delocalized_shell_indices,
        extended_shell_indices=extended_shell_indices,
        unmapped_indices=unmapped_indices,
    )

    return out_smirks


# ================================================================
# Step 2 helpers: reactive center detection
# ================================================================

def _compute_reactive_centers(reaction):
    """Return (reactant_sets, product_sets) of reactive centers (bonding changes) by template. Includes directly covalently attached unmapped neighbors for boundary stability; the covalent shell is built in Step 3."""
    changed_map_numbers = _identify_changed_map_numbers(reaction)
    reactive_centers = ([], [])

    # Reactants
    for i in range(reaction.GetNumReactantTemplates()):
        mol = reaction.GetReactantTemplate(i)
        indices = set()
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() in changed_map_numbers:
                indices.add(atom.GetIdx())
        reactive_centers[0].append(indices)

    # Products
    for i in range(reaction.GetNumProductTemplates()):
        mol = reaction.GetProductTemplate(i)
        indices = set()
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() in changed_map_numbers:
                indices.add(atom.GetIdx())
        reactive_centers[1].append(indices)

    # Add adjacent unmapped atoms on both sides (symmetric)
    for i in range(reaction.GetNumReactantTemplates()):
        mol = reaction.GetReactantTemplate(i)
        for idx in list(reactive_centers[0][i]):
            a = mol.GetAtomWithIdx(idx)
            for n in a.GetNeighbors():
                if n.GetAtomMapNum() == 0:
                    reactive_centers[0][i].add(n.GetIdx())

    for i in range(reaction.GetNumProductTemplates()):
        mol = reaction.GetProductTemplate(i)
        for idx in list(reactive_centers[1][i]):
            a = mol.GetAtomWithIdx(idx)
            for n in a.GetNeighbors():
                if n.GetAtomMapNum() == 0:
                    reactive_centers[1][i].add(n.GetIdx())

    return reactive_centers


def _identify_changed_map_numbers(reaction):
    """Return a set of atom-map numbers whose local environment changes.

    Signature includes:
      - Neighbor identity: (neighbor atomic number, neighbor map number)
      - Tetrahedral stereo feature: map-anchored 4-list encoding parity
      - Alkene stereo feature: map-anchored "diagonals" identity for stereodefined double bonds
      - Adjacent bond features (non-stereo): (neighbor map number, bond type, is_aromatic)

    If an atom's signature differs between sides, its map number is included.
    """

    def _permutation_parity(from_list, to_list):
        """
        Return 0 if the permutation mapping from_list -> to_list is even, 1 if odd.
        Assumes the lists contain unique elements.
        """
        pos = {v: i for i, v in enumerate(from_list)}
        perm = [pos[v] for v in to_list]
        inv = 0
        n = len(perm)
        for i in range(n):
            for j in range(i + 1, n):
                if perm[i] > perm[j]:
                    inv ^= 1
        return inv

    def _tetra_stereo_feature(atom):
        """
        Return a deterministic 4-list encoding tetrahedral parity in a map-anchored neighbor order,
        or None if unspecified or not representable as a 4-neighbor tetra center.
        """
        tag = atom.GetChiralTag()
        if tag == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            return None
        if tag not in (
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
        ):
            return None

        nbr_maps = [nbr.GetAtomMapNum() for nbr in atom.GetNeighbors()]
        if len(nbr_maps) != 4:
            return None

        canon = sorted(nbr_maps)
        rdkit_bit = 0 if tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW else 1
        perm_odd = _permutation_parity(canon, nbr_maps)
        norm_bit = rdkit_bit ^ perm_odd

        if norm_bit == 0:
            return canon
        return [canon[0], canon[1], canon[3], canon[2]]

    def _alkene_diagonals_feature_for_bond(mol, bond):
        """
        If `bond` is a stereodefined double bond, return a canonical "diagonals identity"
        encoded using atom-map numbers:

            (partner_map, ((a,b),(c,d)))

        where ((a,b),(c,d)) are the two diagonal pairs as sorted 2-tuples, ordered
        canonically. Returns None if stereo is unspecified or the bond does not fit
        the 2-substituent-per-end model.
        """
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            return None
        if bond.GetIsAromatic():
            return None

        stereo = bond.GetStereo()
        if stereo not in (Chem.rdchem.BondStereo.STEREOCIS, Chem.rdchem.BondStereo.STEREOTRANS):
            return None

        stereo_atoms = list(bond.GetStereoAtoms())
        if len(stereo_atoms) != 2:
            return None

        b = bond.GetBeginAtom()
        e = bond.GetEndAtom()
        bidx = b.GetIdx()
        eidx = e.GetIdx()

        s0 = mol.GetAtomWithIdx(stereo_atoms[0])
        s1 = mol.GetAtomWithIdx(stereo_atoms[1])

        # Determine which stereo atom is attached to which end.
        def _attached_to_end(s_atom, end_idx, other_end_idx):
            if s_atom.GetIdx() == other_end_idx:
                return False
            return mol.GetBondBetweenAtoms(s_atom.GetIdx(), end_idx) is not None

        b_ref = None
        e_ref = None

        if _attached_to_end(s0, bidx, eidx):
            b_ref = s0
        elif _attached_to_end(s0, eidx, bidx):
            e_ref = s0

        if _attached_to_end(s1, bidx, eidx):
            b_ref = s1
        elif _attached_to_end(s1, eidx, bidx):
            e_ref = s1

        if b_ref is None or e_ref is None:
            return None

        # Collect the two substituents on each end (excluding the partner on the double bond)
        b_subs = [n for n in b.GetNeighbors() if n.GetIdx() != eidx]
        e_subs = [n for n in e.GetNeighbors() if n.GetIdx() != bidx]
        if len(b_subs) != 2 or len(e_subs) != 2:
            return None

        # Identify the "other" substituent at each end
        b_other = b_subs[0] if b_subs[1].GetIdx() == b_ref.GetIdx() else b_subs[1]
        e_other = e_subs[0] if e_subs[1].GetIdx() == e_ref.GetIdx() else e_subs[1]

        # Convert to map numbers (map 0 is allowed and stays as 0)
        mb_ref = b_ref.GetAtomMapNum()
        mb_other = b_other.GetAtomMapNum()
        me_ref = e_ref.GetAtomMapNum()
        me_other = e_other.GetAtomMapNum()

        # Diagonals depend on cis vs trans
        if stereo == Chem.rdchem.BondStereo.STEREOTRANS:
            d1 = tuple(sorted((mb_ref, me_ref)))
            d2 = tuple(sorted((mb_other, me_other)))
        else:
            d1 = tuple(sorted((mb_ref, me_other)))
            d2 = tuple(sorted((mb_other, me_ref)))

        # Canonicalize the set-of-two-sets as an ordered pair
        pair = tuple(sorted((d1, d2)))

        # Encode relative to the bond partner for disambiguation in the atom-level signature
        # If this feature is attached to the begin atom, partner is end's map, and vice versa.
        # Caller will wrap with the appropriate partner map.
        return pair

    def _alkene_features_for_atom(atom):
        """
        Return a frozenset of alkene stereo features for stereodefined double bonds incident to `atom`.
        Each entry is:
            (partner_map, ((a,b),(c,d)))
        where ((a,b),(c,d)) is the canonical diagonals identity.
        Returns None if no stereodefined double bonds incident, or if stereo is absent.
        """
        feats = set()
        mol = atom.GetOwningMol()
        aidx = atom.GetIdx()

        for bond in atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue
            if bond.GetIsAromatic():
                continue

            other = bond.GetOtherAtom(atom)
            partner_map = other.GetAtomMapNum()

            diag_pair = _alkene_diagonals_feature_for_bond(mol, bond)
            if diag_pair is None:
                continue

            feats.add((partner_map, diag_pair))

        if not feats:
            return None
        return frozenset(feats)

    def _atom_signature(atom):
        neighbors = set()
        bond_features = set()

        tetra = _tetra_stereo_feature(atom)
        alkene = _alkene_features_for_atom(atom)

        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            nbr_map = nbr.GetAtomMapNum()

            neighbors.add((nbr.GetAtomicNum(), nbr_map))

            btype = int(bond.GetBondType())
            aromatic = bond.GetIsAromatic()

            # Important: do NOT include RDKit bond stereo enum here, it is not map-anchored.
            bond_features.add((nbr_map, btype, aromatic))

        return {
            "neighbors": neighbors,
            "tetra": tuple(tetra) if tetra is not None else None,
            "alkene": tuple(sorted(alkene)) if alkene is not None else None,
            "bond_features": bond_features,
        }

    def _signatures_for(mol):
        out = {}
        for atom in mol.GetAtoms():
            amap = atom.GetAtomMapNum()
            if amap > 0:
                out[amap] = _atom_signature(atom)
        return out

    reactant_sigs, product_sigs = {}, {}
    for i in range(reaction.GetNumReactantTemplates()):
        rtpl = reaction.GetReactantTemplate(i)
        reactant_sigs.update(_signatures_for(rtpl))
    for i in range(reaction.GetNumProductTemplates()):
        ptpl = reaction.GetProductTemplate(i)
        product_sigs.update(_signatures_for(ptpl))

    changed = set()
    for amap in set(reactant_sigs) | set(product_sigs):
        if reactant_sigs.get(amap) != product_sigs.get(amap):
            changed.add(amap)
    return changed


# ================================================================
# Step 3 helper: covalent (sigma) shell
# ================================================================

def _compute_covalent_shells(reaction, reactive_centers):
    """Return (reactant_sets, product_sets) of covalent-shell atoms (σ) by template."""
    covalent_shell_indices = ([], [])
    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        covalent_shell_indices[0].append(
            _collect_covalent_shell(reactant, reactive_centers[0][i])
        )
    for i in range(reaction.GetNumProductTemplates()):
        product = reaction.GetProductTemplate(i)
        covalent_shell_indices[1].append(
            _collect_covalent_shell(product, reactive_centers[1][i])
        )
    return covalent_shell_indices


def _collect_covalent_shell(molecule, reactive_center_indices):
    """Return covalent-shell atoms (σ): those one bond from any reactive center in the given molecule."""
    covalent_shell_indices = set()
    for atom in molecule.GetAtoms():
        atom_idx = atom.GetIdx()
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in reactive_center_indices:
                covalent_shell_indices.add(atom_idx)
    return covalent_shell_indices


# ================================================================
# Step 4 helpers: delocalized (pi) growth
# ================================================================

def _grow_delocalized_shell(reaction, covalent_shell_indices):
    """Build the delocalized shell (π): iterative growth from the σ shell across conjugation/unsaturation."""

    def _state_key(atom_sets):
        return (
            tuple(frozenset(s) for s in atom_sets[0]),
            tuple(frozenset(s) for s in atom_sets[1]),
        )

    def _expand_molecule(mol, covalent_shell_atoms, mapped_ids):
        delocalized_atoms = set()
        for idx in covalent_shell_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if _is_conjugation_capable(atom):
                delocalized_atoms.add(idx)

                amap = atom.GetAtomMapNum()
                if amap > 0:
                    mapped_ids.add(amap)

                for nbr in atom.GetNeighbors():
                    if _is_conjugation_capable(nbr):
                        delocalized_atoms.add(nbr.GetIdx())

        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() in mapped_ids:
                delocalized_atoms.add(atom.GetIdx())

        return delocalized_atoms

    prev = (
        [set(s) for s in covalent_shell_indices[0]],
        [set(s) for s in covalent_shell_indices[1]],
    )

    while True:
        mapped_ids = set()
        new_reactant_sets, new_product_sets = [], []

        for i in range(reaction.GetNumReactantTemplates()):
            mol = reaction.GetReactantTemplate(i)
            new_reactant_sets.append(_expand_molecule(mol, prev[0][i], mapped_ids))

        for i in range(reaction.GetNumProductTemplates()):
            mol = reaction.GetProductTemplate(i)
            new_product_sets.append(_expand_molecule(mol, prev[1][i], mapped_ids))

        for i in range(reaction.GetNumReactantTemplates()):
            mol = reaction.GetReactantTemplate(i)
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() in mapped_ids:
                    new_reactant_sets[i].add(atom.GetIdx())

        for i in range(reaction.GetNumProductTemplates()):
            mol = reaction.GetProductTemplate(i)
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() in mapped_ids:
                    new_product_sets[i].add(atom.GetIdx())

        if _state_key((new_reactant_sets, new_product_sets)) == _state_key(prev):
            return (new_reactant_sets, new_product_sets)

        prev = (new_reactant_sets, new_product_sets)


def _is_conjugation_capable(atom):
    """Return True if any bond is not single (proxy for conjugation/aromaticity)."""
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return True
    return False


# ================================================================
# Step 5 helper: extended shell (sigma neighbors of delocalized shell)
# ================================================================

def _add_extended_shell(reaction, covalent_shell_indices, delocalized_shell_indices):
    """Return σ-neighbors of π shell, excluding existing covalent/delocalized members."""
    extended = ([], [])

    for i in range(reaction.GetNumReactantTemplates()):
        mol = reaction.GetReactantTemplate(i)
        cov = set(covalent_shell_indices[0][i])
        dloc = set(delocalized_shell_indices[0][i])
        add = set()
        for idx in dloc:
            a = mol.GetAtomWithIdx(idx)
            for nbr in a.GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx not in cov and nidx not in dloc:
                    add.add(nidx)
        extended[0].append(add)

    for i in range(reaction.GetNumProductTemplates()):
        mol = reaction.GetProductTemplate(i)
        cov = set(covalent_shell_indices[1][i])
        dloc = set(delocalized_shell_indices[1][i])
        add = set()
        for idx in dloc:
            a = mol.GetAtomWithIdx(idx)
            for nbr in a.GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx not in cov and nidx not in dloc:
                    add.add(nidx)
        extended[1].append(add)

    return extended


# ================================================================
# Step 6 helper: unmapped atoms
# ================================================================

def _compute_unmapped_sets(reaction):
    """Return (reactant_sets, product_sets) of unmapped atom indices by template."""
    unmapped_indices = ([], [])
    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        unmapped_indices[0].append(_collect_unmapped_atoms(reactant))
    for i in range(reaction.GetNumProductTemplates()):
        product = reaction.GetProductTemplate(i)
        unmapped_indices[1].append(_collect_unmapped_atoms(product))
    return unmapped_indices


def _collect_unmapped_atoms(molecule):
    """Return set of unmapped atom indices for a molecule (element-agnostic)."""
    unmapped_indices = set()
    for atom in molecule.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            unmapped_indices.add(atom.GetIdx())
    return unmapped_indices


# ================================================================
# Step 7: Stereochemistry utilities
# ================================================================

def _strip_stereochemistry(molecule: Chem.Mol) -> Chem.Mol:
    for atom in molecule.GetAtoms():
        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
    for bond in molecule.GetBonds():
        bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)
    return molecule


# ================================================================
# Step 8: compile operator
# ================================================================

def _compile_operator(
    *,
    reaction,
    reactive_centers,
    covalent_shell_indices,
    delocalized_shell_indices,
    extended_shell_indices,
    unmapped_indices,
):
    keep_atom_indices = ([], [])
    for i in range(len(reactive_centers[0])):
        keep_atom_indices[0].append(set(reactive_centers[0][i]))
    for i in range(len(reactive_centers[1])):
        keep_atom_indices[1].append(set(reactive_centers[1][i]))

    if covalent_shell_indices is not None:
        for i in range(len(covalent_shell_indices[0])):
            keep_atom_indices[0][i].update(covalent_shell_indices[0][i])
        for i in range(len(covalent_shell_indices[1])):
            keep_atom_indices[1][i].update(covalent_shell_indices[1][i])

    if delocalized_shell_indices is not None:
        for i in range(len(delocalized_shell_indices[0])):
            keep_atom_indices[0][i].update(delocalized_shell_indices[0][i])
        for i in range(len(delocalized_shell_indices[1])):
            keep_atom_indices[1][i].update(delocalized_shell_indices[1][i])

    if extended_shell_indices is not None:
        for i in range(len(extended_shell_indices[0])):
            keep_atom_indices[0][i].update(extended_shell_indices[0][i])
        for i in range(len(extended_shell_indices[1])):
            keep_atom_indices[1][i].update(extended_shell_indices[1][i])

    if unmapped_indices is not None:
        for i in range(len(unmapped_indices[0])):
            keep_atom_indices[0][i].update(unmapped_indices[0][i])
        for i in range(len(unmapped_indices[1])):
            keep_atom_indices[1][i].update(unmapped_indices[1][i])

    remove_bond_indices = ([], [])
    remove_atom_indices = ([], [])

    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        bond_indices_set = set()
        atom_indices_set = set()
        for atom in reactant.GetAtoms():
            if atom.GetIdx() not in keep_atom_indices[0][i]:
                atom_indices_set.add(atom.GetIdx())
                for bond in atom.GetBonds():
                    bond_indices_set.add(bond.GetIdx())
        remove_bond_indices[0].append(bond_indices_set)
        remove_atom_indices[0].append(atom_indices_set)

    for i in range(reaction.GetNumProductTemplates()):
        product = reaction.GetProductTemplate(i)
        bond_indices_set = set()
        atom_indices_set = set()
        for atom in product.GetAtoms():
            if atom.GetIdx() not in keep_atom_indices[1][i]:
                atom_indices_set.add(atom.GetIdx())
                for bond in atom.GetBonds():
                    bond_indices_set.add(bond.GetIdx())
        remove_bond_indices[1].append(bond_indices_set)
        remove_atom_indices[1].append(atom_indices_set)

    new_reaction = rdChemReactions.ChemicalReaction()

    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        editable_reactant = Chem.EditableMol(reactant)
        for bond_idx in remove_bond_indices[0][i]:
            bond = reactant.GetBondWithIdx(bond_idx)
            editable_reactant.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for atom_idx in sorted(remove_atom_indices[0][i], reverse=True):
            editable_reactant.RemoveAtom(atom_idx)
        r_mol = editable_reactant.GetMol()
        new_reaction.AddReactantTemplate(r_mol)

    for i in range(reaction.GetNumProductTemplates()):
        product = reaction.GetProductTemplate(i)
        editable_product = Chem.EditableMol(product)
        for bond_idx in remove_bond_indices[1][i]:
            bond = product.GetBondWithIdx(bond_idx)
            editable_product.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for atom_idx in sorted(remove_atom_indices[1][i], reverse=True):
            editable_product.RemoveAtom(atom_idx)
        p_mol = editable_product.GetMol()
        new_reaction.AddProductTemplate(p_mol)

    _unset_hydrogen_map_on_noncarbon(new_reaction)

    return rdChemReactions.ReactionToSmarts(new_reaction)


def extract_operator_by_abstraction(smirks: str, abstraction: str, matched: bool = False):
    """Dispatch into `operator_extractor` using a discrete EVODEX abstraction level."""
    level = abstraction.upper()

    if level == "A":
        b_operator = operator_extractor(
            smirks,
            include_stereochemistry=True,
            include_sigma=False,
            include_pi=False,
            include_extended=False,
            include_unmapped=False,
        )
        return _abstract_operator_atoms_to_wildcards(b_operator)

    if level == "B":
        return operator_extractor(
            smirks,
            include_stereochemistry=True,
            include_sigma=False,
            include_pi=False,
            include_extended=False,
            include_unmapped=False,
        )
    elif level == "C":
        return operator_extractor(
            smirks,
            include_stereochemistry=True,
            include_sigma=True,
            include_pi=False,
            include_extended=False,
            include_unmapped=False,
        )
    elif level == "D":
        return operator_extractor(
            smirks,
            include_stereochemistry=True,
            include_sigma=True,
            include_pi=True,
            include_extended=False,
            include_unmapped=False,
        )
    elif level == "E":
        return operator_extractor(
            smirks,
            include_stereochemistry=True,
            include_sigma=True,
            include_pi=True,
            include_extended=True,
            include_unmapped=False,
        )


# ================================================================
# Helper: Abstract atom identity to wildcards for level A
# ================================================================

def _abstract_operator_atoms_to_wildcards(operator_smirks: str) -> str:
    """Return a SMIRKS where all atoms have wildcard identity but keep map numbers."""
    s = operator_smirks
    out = []
    i = 0
    n = len(s)

    while i < n:
        ch = s[i]

        if ch == "[":
            j = s.find("]", i + 1)
            if j == -1:
                out.append(s[i:])
                break
            interior = s[i + 1:j]
            m = re.search(r":(\d+)", interior)
            if m:
                out.append(f"[*:{m.group(1)}]")
            else:
                out.append("*")
            i = j + 1

        elif ch.isalpha():
            if ch.isupper() and i + 1 < n and s[i + 1].islower():
                i += 2
            else:
                i += 1
            out.append("*")

        else:
            out.append(ch)
            i += 1

    return "".join(out)


# ================================================================
# Helper: Remove map numbers from hydrogen atoms attached to non-carbon
# ================================================================

def _unset_hydrogen_map_on_noncarbon(reaction):
    """
    Remove map numbers from hydrogen-like atoms (Z=1 or 85) attached to non-carbon atoms.
    """
    maps_to_clear = set()

    def _collect_maps_to_clear(mol):
        for atom in mol.GetAtoms():
            amap = atom.GetAtomMapNum()
            if amap == 0:
                continue

            z = atom.GetAtomicNum()
            if z not in (1, 85):
                continue

            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 6:
                    maps_to_clear.add(amap)
                    break

    for i in range(reaction.GetNumReactantTemplates()):
        _collect_maps_to_clear(reaction.GetReactantTemplate(i))

    for i in range(reaction.GetNumProductTemplates()):
        _collect_maps_to_clear(reaction.GetProductTemplate(i))

    if not maps_to_clear:
        return reaction

    def _clear_maps(mol):
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() in maps_to_clear:
                atom.SetAtomMapNum(0)

    for i in range(reaction.GetNumReactantTemplates()):
        _clear_maps(reaction.GetReactantTemplate(i))

    for i in range(reaction.GetNumProductTemplates()):
        _clear_maps(reaction.GetProductTemplate(i))

    return reaction


if __name__ == "__main__":
    examples = {
        "alcohol_methylation": "[H][O:2][C:3]([H:4])([H:5])[C:6]([H:7])([H:8])[H:9]>>[C]([H])([H])([H])[O:2][C:3]([H:4])([H:5])[C:6]([H:7])([H:8])[H:9]",
        "phenol_methylation": "[H][O:2][c:3]1[c:4]([H:5])[c:6]([H:7])[c:8]([H:9])[c:10]([H:11])[c:12]1[H:13]>>[C]([H])([H])([H])[O:2][c:3]1[c:4]([H:5])[c:6]([H:7])[c:8]([H:9])[c:10]([H:11])[c:12]1[H:13]",
        "cis_trans_isomerization": "[C:1]([H:4])([H:5])([H:6])/[C:2]([H:7])=[C:3]([H:8])/[C:4]([H:9])([H:10])([H:11])>>[C:1]([H:4])([H:5])([H:6])\\[C:2]([H:7])=[C:3]([H:8])/[C:4]([H:9])([H:10])([H:11])",
        "nucleophilic_aromatic_substitution": "[C:1]([c:2]1[c:3]([H:12])[c:4]([H:13])[c:5]([F])[c:7]([H:14])[c:8]1[H:15])([H:9])([H:10])[H:11]>>[C:1]([c:2]1[c:3]([H:12])[c:4]([H:13])[c:5]([Cl])[c:7]([H:14])[c:8]1[H:15])([H:9])([H:10])[H:11]",

        # --------------------------
        # Tetrahedral stereo tests
        # --------------------------
        "tetra_inversion": "[C@:1]([F:2])([Cl:3])([Br:4])[I:5]>>[C@@:1]([F:2])([Cl:3])([Br:4])[I:5]",
        "tetra_at_at_reordered": "[C@:1]([F:2])([Cl:3])([Br:4])[I:5]>>[C@:1]([Cl:3])([F:2])([Br:4])[I:5]",
        "tetra_same_config_reordered": "[C@:1]([F:2])([Cl:3])([Br:4])[I:5]>>[C@@:1]([Cl:3])([F:2])([Br:4])[I:5]",

        # --------------------------
        # Alkene E/Z stereo tests
        # --------------------------

        # True E/Z change (cis <-> trans) with explicit H substituents.
        # This should be detected as a change at the alkene carbons.
        "ez_isomerization": "[F:3]/[C:1]([H:4])=[C:2]([H:6])/[Cl:5]>>[F:3]/[C:1]([H:4])=[C:2]([H:6])\\[Cl:5]",

        # Same absolute configuration, both slashes flipped (should be unchanged).
        "ez_same_config_flip_both": "[F:3]/[C:1]([H:4])=[C:2]([H:6])\\[Cl:5]>>[F:3]\\[C:1]([H:4])=[C:2]([H:6])/[Cl:5]",

        # Same absolute configuration, reorder branch presentation around the alkene end (should be unchanged).
        "ez_same_config_reordered_text": "[F:3]/[C:1]([H:4])=[C:2]([H:6])\\[Cl:5]>>[C:1]([H:4])(/[F:3])=[C:2]([H:6])\\[Cl:5]",
    }

    for name, smirks in examples.items():
        print(f"\nExample: {name}")
        print(f"  Full SMIRKS: {smirks}")
        rxn = rdChemReactions.ReactionFromSmarts(smirks)
        rxn.Initialize()
        print(f"  Changed map numbers: {sorted(_identify_changed_map_numbers(rxn))}")
        op_A = extract_operator_by_abstraction(smirks, "A")
        print(f"  A abstraction: {op_A}")