from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

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
    # Build the reaction used for compilation, stripping stereochemistry if requested
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


    # Step 8: Compile operator: assemble keep sets, prune others, and emit SMIRKS
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
    """Return a set of atom-map numbers whose local bonding environment changes.

    For each mapped atom on the reactant and product sides we build a *signature*
    that captures any feature that could signal a bonding change:
      - Neighbor identity: (neighbor atomic number, neighbor map number)
      - Atom chirality: the atom's RDKit `ChiralTag` name
      - Adjacent bond features for each incident bond:
          (neighbor map number, bond type, bond stereo, is_aromatic)

    If an atom's signature differs between sides, its map number is included in
    the result set.
    """

    def _atom_signature(atom):
        """Construct a comparison-ready signature for a single atom."""
        neighbors = set()
        bond_features = set()
        chirality = atom.GetChiralTag().name

        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            nbr_map = nbr.GetAtomMapNum()

            # Neighbor identity
            neighbors.add((nbr.GetAtomicNum(), nbr_map))

            # Bond attributes that affect E/Z and conjugation context
            btype = int(bond.GetBondType())
            stereo_enum = bond.GetStereo()
            stereo = getattr(stereo_enum, 'name', int(stereo_enum))
            aromatic = bond.GetIsAromatic()

            # Key by neighbor map number to keep it deterministic
            bond_features.add((nbr_map, btype, stereo, aromatic))

        return {
            'neighbors': neighbors,
            'chirality': chirality,
            'bond_features': bond_features,
        }

    def _signatures_for(mol):
        """Map: atom map number -> atom signature (only for mapped atoms)."""
        out = {}
        for atom in mol.GetAtoms():
            amap = atom.GetAtomMapNum()
            if amap > 0:
                out[amap] = _atom_signature(atom)
        return out

    # Build signatures for all templates on each side
    reactant_sigs, product_sigs = {}, {}
    for i in range(reaction.GetNumReactantTemplates()):
        rtpl = reaction.GetReactantTemplate(i)
        reactant_sigs.update(_signatures_for(rtpl))
    for i in range(reaction.GetNumProductTemplates()):
        ptpl = reaction.GetProductTemplate(i)
        product_sigs.update(_signatures_for(ptpl))

    # Compare signatures; any difference means a bonding change
    changed = set()
    for amap in set(reactant_sigs) | set(product_sigs):
        if reactant_sigs.get(amap) != product_sigs.get(amap):
            changed.add(amap)
    return changed


# ================================================================
# Step 3 helper: covalent (sigma) shell
# ================================================================

def _compute_covalent_shells(reaction, reactive_centers):
    """Return (reactant_sets, product_sets) of covalent-shell atoms (σ) by template.
    Thin wrapper that aggregates per-molecule results from `_collect_covalent_shell`.
    """
    covalent_shell_indices = ([], [])
    # Reactants
    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        covalent_shell_indices[0].append(
            _collect_covalent_shell(reactant, reactive_centers[0][i])
        )
    # Products
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
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in reactive_center_indices:
                covalent_shell_indices.add(atom_idx)
    return covalent_shell_indices


# ================================================================
# Step 4 helpers: delocalized (pi) growth
# ================================================================

def _grow_delocalized_shell(reaction, covalent_shell_indices):
    """Build the delocalized shell (π): iterative growth from the σ shell across conjugation/unsaturation; mirrored via atom-map numbers; excludes inductive and sp3-only participation by design of _is_conjugation_capable."""

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
    """Return True if any bond is not single (proxy for conjugation/aromaticity). Intentionally ignores purely inductive pathways and sp3-only participation."""
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return True
    return False

# ================================================================
# Step 5 helper: extended shell (sigma neighbors of delocalized shell)
# ================================================================

def _add_extended_shell(reaction, covalent_shell_indices, delocalized_shell_indices):
    """Return (reactant_sets, product_sets) of atoms that are σ-neighbors of any delocalized-shell atom,
    excluding atoms already in the covalent or delocalized shells.
    """
    extended = ([], [])

    # Reactants
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

    # Products
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
    # Reactants
    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        unmapped_indices[0].append(_collect_unmapped_atoms(reactant))
    # Products
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
# Assemble keep sets, prune atoms and bonds not kept, and emit SMIRKS
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
    """
    Assemble the operator by selecting which atoms to keep, pruning the rest, and
    emitting a compact SMIRKS.

    Inputs (summarized)
    - reactive_centers: atoms whose bonding changed (per template, reactants/products)
    - covalent_shell_indices: atoms one σ bond from the centers (included if not None)
    - delocalized_shell_indices: conjugation-driven shell grown from σ (included if not None)
    - extended_shell_indices: σ-neighbors of the delocalized shell (included if not None)
    - unmapped_indices: atoms with no map numbers (included if not None)

    Shells are included when their corresponding per-template sets are not None.

    Output
    - SMIRKS string of the pruned reaction.
    """

    # 8.1) Initialize per-template keep-sets with reactive centers (always kept)
    keep_atom_indices = ([], [])
    for i in range(len(reactive_centers[0])):
        keep_atom_indices[0].append(set(reactive_centers[0][i]))
    for i in range(len(reactive_centers[1])):
        keep_atom_indices[1].append(set(reactive_centers[1][i]))

    # 8.2) Include the covalent (σ) shell if present
    if covalent_shell_indices is not None:
        for i in range(len(covalent_shell_indices[0])):
            keep_atom_indices[0][i].update(covalent_shell_indices[0][i])
        for i in range(len(covalent_shell_indices[1])):
            keep_atom_indices[1][i].update(covalent_shell_indices[1][i])

    # 8.3) Include the delocalized (π) shell if present
    if delocalized_shell_indices is not None:
        for i in range(len(delocalized_shell_indices[0])):
            keep_atom_indices[0][i].update(delocalized_shell_indices[0][i])
        for i in range(len(delocalized_shell_indices[1])):
            keep_atom_indices[1][i].update(delocalized_shell_indices[1][i])

    # 8.4) Include the extended shell (σ neighbors of π shell) if present
    if extended_shell_indices is not None:
        for i in range(len(extended_shell_indices[0])):
            keep_atom_indices[0][i].update(extended_shell_indices[0][i])
        for i in range(len(extended_shell_indices[1])):
            keep_atom_indices[1][i].update(extended_shell_indices[1][i])

    # 8.5) Include unmapped atoms (appear on only one side) if present
    if unmapped_indices is not None:
        for i in range(len(unmapped_indices[0])):
            keep_atom_indices[0][i].update(unmapped_indices[0][i])
        for i in range(len(unmapped_indices[1])):
            keep_atom_indices[1][i].update(unmapped_indices[1][i])

    # 8.6) Compute removal sets per template
    #       Any atom not in the keep-set is marked for deletion and all incident bonds
    #       are marked for removal.
    remove_bond_indices = ([], [])
    remove_atom_indices = ([], [])

    # 8.6a) Reactant-side removals
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

    # 8.6b) Product-side removals
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

    # 8.7) Construct a new reaction from pruned templates
    new_reaction = rdChemReactions.ChemicalReaction()

    # 8.7a) Rebuild reactant templates with unwanted atoms/bonds pruned
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

    # 8.7b) Rebuild product templates with unwanted atoms/bonds pruned
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

    return rdChemReactions.ReactionToSmarts(new_reaction)


if __name__ == "__main__":
    # Alcohol Methylation
    # smirks = "[H][O:2][C:3]([H:4])([H:5])[C:6]([H:7])([H:8])[H:9]>>[C]([H])([H])([H])[O:2][C:3]([H:4])([H:5])[C:6]([H:7])([H:8])[H:9]"
    
    # Phenol Methylation
    # smirks = '[H][O:2][c:3]1[c:4]([H:5])[c:6]([H:7])[c:8]([H:9])[c:10]([H:11])[c:12]1[H:13]>>[C]([H])([H])([H])[O:2][c:3]1[c:4]([H:5])[c:6]([H:7])[c:8]([H:9])[c:10]([H:11])[c:12]1[H:13]'
    
    # cis/trans Isomerization
    # smirks = "[C:1]([H:4])([H:5])([H:6])/[C:2]([H:7])=[C:3]([H:8])/[C:4]([H:9])([H:10])([H:11])>>[C:1]([H:4])([H:5])([H:6])\\[C:2]([H:7])=[C:3]([H:8])/[C:4]([H:9])([H:10])([H:11])"

    # Nucleophillic Aromatic Substitution
    smirks = "[C:1]([c:2]1[c:3]([H:12])[c:4]([H:13])[c:5]([F])[c:7]([H:14])[c:8]1[H:15])([H:9])([H:10])[H:11]>>[C:1]([c:2]1[c:3]([H:12])[c:4]([H:13])[c:5]([Cl])[c:7]([H:14])[c:8]1[H:15])([H:9])([H:10])[H:11]"
    # smirks = "[CH3:1][C:2]([CH3:3])=[CH:4][CH2:5][O:24][P:22](=[O:21])([OH:23])[O:25][P:26](=[O:27])([OH:28])[OH:29].[cH:6]1[cH:7][cH:8][cH:9][c:10]2[nH:11][cH:12][c:13]([CH2:14][C@H:15]([NH2:16])[C:17](=[O:18])[OH:19])[c:20]12>>[CH3:1][C:2]([CH3:3])=[CH:4][CH2:5][c:6]1[cH:7][cH:8][cH:9][c:10]2[nH:11][cH:12][c:13]([CH2:14][C@H:15]([NH2:16])[C:17](=[O:18])[OH:19])[c:20]12.[O:21]=[P:22]([OH:23])([OH:24])[O:25][P:26](=[O:27])([OH:28])[OH:29]"

    result = operator_extractor(
        smirks,
        include_stereochemistry=True,
        include_sigma=True,
        include_pi=True,
        include_extended=True,
        include_unmapped=True,
    )
    print(result)