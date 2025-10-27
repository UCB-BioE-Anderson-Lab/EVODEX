from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Function to remove stereochemistry from a molecule
def _remove_stereochemistry(molecule):
    for atom in molecule.GetAtoms():
        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
    for bond in molecule.GetBonds():
        bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)
    return molecule

# Compute center atom indices on both sides from changed map numbers and include adjacent unmapped atoms
# Returns: (reactant_center_sets, product_center_sets)
def _compute_center_atom_indices(reaction):
    center_atom_maps = _identify_changed_map_numbers(reaction)
    center_atom_indices = ([], [])

    # Reactants
    for i in range(reaction.GetNumReactantTemplates()):
        mol = reaction.GetReactantTemplate(i)
        indices = set()
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() in center_atom_maps:
                indices.add(atom.GetIdx())
        center_atom_indices[0].append(indices)

    # Products
    for i in range(reaction.GetNumProductTemplates()):
        mol = reaction.GetProductTemplate(i)
        indices = set()
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() in center_atom_maps:
                indices.add(atom.GetIdx())
        center_atom_indices[1].append(indices)

    # Add adjacent unmapped atoms on both sides (symmetric)
    for i in range(reaction.GetNumReactantTemplates()):
        mol = reaction.GetReactantTemplate(i)
        for idx in list(center_atom_indices[0][i]):
            a = mol.GetAtomWithIdx(idx)
            for n in a.GetNeighbors():
                if n.GetAtomMapNum() == 0:
                    center_atom_indices[0][i].add(n.GetIdx())

    for i in range(reaction.GetNumProductTemplates()):
        mol = reaction.GetProductTemplate(i)
        for idx in list(center_atom_indices[1][i]):
            a = mol.GetAtomWithIdx(idx)
            for n in a.GetNeighbors():
                if n.GetAtomMapNum() == 0:
                    center_atom_indices[1][i].add(n.GetIdx())

    return center_atom_indices

# Identify atom-map numbers whose local environments differ between reactants and products
# Returns: set[int] of changed atom-map numbers
def _identify_changed_map_numbers(reaction):
    def get_atom_signature(atom):
        sig = {'neighbors': set(), 'chirality': atom.GetChiralTag().name}
        for bond in atom.GetBonds():
            n = bond.GetOtherAtom(atom)
            sig['neighbors'].add((n.GetAtomicNum(), n.GetAtomMapNum()))
        return sig

    def build_signature_dict(mol):
        out = {}
        for a in mol.GetAtoms():
            amap = a.GetAtomMapNum()
            if amap > 0:
                out[amap] = get_atom_signature(a)
        return out

    reactant_sigs, product_sigs = {}, {}
    for i in range(reaction.GetNumReactantTemplates()):
        reactant_sigs.update(build_signature_dict(reaction.GetReactantTemplate(i)))
    for i in range(reaction.GetNumProductTemplates()):
        product_sigs.update(build_signature_dict(reaction.GetProductTemplate(i)))

    changed = set()
    for amap in set(reactant_sigs) | set(product_sigs):
        if reactant_sigs.get(amap) != product_sigs.get(amap):
            changed.add(amap)
    return changed

# Function to process a molecule and identify center atom indices
def _extract_center_atom_indices(molecule, center_atom_maps):
    center_atom_indices = set()
    for atom in molecule.GetAtoms():
        atom_map = atom.GetAtomMapNum()
        if atom_map in center_atom_maps:
            center_atom_indices.add(atom.GetIdx())
    return center_atom_indices

# Function to process a molecule and identify static hydrogen atom indices
def _process_static_hydrogens(molecule, center_atom_indices):
    static_hydrogen_indices = set()
    for atom in molecule.GetAtoms():
        # Skip if not a hydrogen
        if atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 85: 
            continue
        # Skip if  unmapped
        if atom.GetAtomMapNum() == 0: 
            continue
        # Skip if its a reactive center
        if atom.GetIdx() not in center_atom_indices:
            static_hydrogen_indices.add(atom.GetIdx())
    return static_hydrogen_indices

# Function to process a molecule and identify sigma atom indices
# Identifies all atoms adjacent to reactive centers (radius = 1)
def _process_sigma_molecule(molecule, center_atom_indices):
    sigma_indices = set()
    for atom in molecule.GetAtoms():
        atom_idx = atom.GetIdx()
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in center_atom_indices:
                sigma_indices.add(atom_idx)
    return sigma_indices

# Helper function to check if an atom's hybridization state supports conjugation
# Implemented as identifying if an atom is not saturated.  If any of the bonds are not single bonds,
# then there is either a double, triple, or aromatic bond.
def _is_pi_capable(atom):
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return True
    return False

 # Iterative breadth-first expansion through the π system
def _grow_pi_shell(reaction, sigma_atom_indices):
    """Expand the π shell from the σ shell across all reactants and products."""

    def _state_key(atom_sets):
        # Order-independent snapshot used to detect when nothing new is added
        return (
            tuple(frozenset(s) for s in atom_sets[0]),
            tuple(frozenset(s) for s in atom_sets[1]),
        )

    def _expand_molecule(mol, sigma_atoms, mapped_ids):
        # From σ atoms in one molecule, keep π-capable seeds, add π-capable neighbors,
        # and record map numbers to mirror counterparts across molecules.
        pi_atoms = set()
        for idx in sigma_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if _is_pi_capable(atom):
                pi_atoms.add(idx)

                amap = atom.GetAtomMapNum()
                if amap > 0:
                    mapped_ids.add(amap)

                for nbr in atom.GetNeighbors():
                    if _is_pi_capable(nbr):
                        pi_atoms.add(nbr.GetIdx())

        # Add any atoms whose map number matches those seen this round
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() in mapped_ids:
                pi_atoms.add(atom.GetIdx())

        return pi_atoms

    # Start from the σ shell (reactants, products)
    prev = (
        [set(s) for s in sigma_atom_indices[0]],
        [set(s) for s in sigma_atom_indices[1]],
    )

    while True:
        mapped_ids = set()  # atom-map numbers used to synchronize counterparts
        new_reactant_sets, new_product_sets = [], []

        # Reactants
        for i in range(reaction.GetNumReactantTemplates()):
            mol = reaction.GetReactantTemplate(i)
            new_reactant_sets.append(_expand_molecule(mol, prev[0][i], mapped_ids))

        # Products
        for i in range(reaction.GetNumProductTemplates()):
            mol = reaction.GetProductTemplate(i)
            new_product_sets.append(_expand_molecule(mol, prev[1][i], mapped_ids))

        # Mirror mapped counterparts across all molecules
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

        # Stop when the sets no longer change
        if _state_key((new_reactant_sets, new_product_sets)) == _state_key(prev):
            return (new_reactant_sets, new_product_sets)

        prev = (new_reactant_sets, new_product_sets)

# Function to process a molecule and identify unmapped hydrogen and heavy atom indices
def _process_unmapped_atoms(molecule):
    unmapped_hydrogen_indices = set()
    unmapped_heavy_indices = set()
    for atom in molecule.GetAtoms():
        atom_map = atom.GetAtomMapNum()
        if atom_map != 0:
            continue
        atomic_number = atom.GetAtomicNum()
        if atomic_number == 1 or atomic_number == 85:
            unmapped_hydrogen_indices.add(atom.GetIdx())
        else:
            unmapped_heavy_indices.add(atom.GetIdx())
    return unmapped_hydrogen_indices, unmapped_heavy_indices

# Main function to process the reaction
def extract_operator(smirks: str, include_stereochemistry: bool = False, include_sigma: bool = True, include_pi: bool = True, include_unmapped_hydrogens: bool = True, include_unmapped_heavy_atoms: bool = True, include_static_hydrogens: bool = False):

    # ==================================================================
    #                   REACTION STEREOCHEMISTRY AND PREP
    # ==================================================================

    # Parse the reaction SMARTS
    reaction = rdChemReactions.ReactionFromSmarts(smirks)

    # Remove stereochemistry if requested
    if not include_stereochemistry:
        # Create a new reaction object
        new_reaction = rdChemReactions.ChemicalReaction()

        # Process reactants
        for i in range(reaction.GetNumReactantTemplates()):
            reactant = reaction.GetReactantTemplate(i)
            reactant = _remove_stereochemistry(reactant)
            new_reaction.AddReactantTemplate(reactant)

        # Process products
        for i in range(reaction.GetNumProductTemplates()):
            product = reaction.GetProductTemplate(i)
            product = _remove_stereochemistry(product)
            new_reaction.AddProductTemplate(product)

        smirks = rdChemReactions.ReactionToSmarts(new_reaction)
        reaction = rdChemReactions.ReactionFromSmarts(smirks)

    # Do some cleanup and internal calculation
    reaction.Initialize()

    # ==================================================================
    #                       LABEL AND COLLECT ATOMS SETS
    # ==================================================================

    # ---------------------- POPULATE CENTER ATOMS ---------------------
    center_atom_indices = _compute_center_atom_indices(reaction)

    # ---------------------- POPULATE SIGMA-BONDED ATOMS ---------------------
    # Create a tuple of lists of sets called sigma_atom_indices
    sigma_atom_indices = ([], [])

    # Iterate through reactants and add to the first list in the tuple
    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        sigma_atom_indices[0].append(_process_sigma_molecule(reactant, center_atom_indices[0][i]))

    # Iterate through products and add to the second list in the tuple
    for i in range(reaction.GetNumProductTemplates()):
        product = reaction.GetProductTemplate(i)
        sigma_atom_indices[1].append(_process_sigma_molecule(product, center_atom_indices[1][i]))

    # ---------------------- POPULATE PI-BONDED ATOMS ---------------------
    # Grow the pi shell starting from the sigma shell
    # The reactive center too may be a Pi bond, but those are already added to the final retained
    # atoms for all settings, and thus they do not need to be traversed again.
    pi_atom_indices = _grow_pi_shell(reaction, sigma_atom_indices)

    # ---------------------- POPULATE UNMAPPED ATOMS ---------------------
    # Initialize unmapped_hydrogen_indices and unmapped_heavy_indices with the same structure as pi_atom_indices
    unmapped_hydrogen_indices = ([], [])
    unmapped_heavy_indices = ([], [])

    # Iterate through reactants and add to the first list in the tuples
    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        hydrogen_indices, heavy_indices = _process_unmapped_atoms(reactant)
        unmapped_hydrogen_indices[0].append(hydrogen_indices)
        unmapped_heavy_indices[0].append(heavy_indices)

    # Iterate through products and add to the second list in the tuples
    for i in range(reaction.GetNumProductTemplates()):
        product = reaction.GetProductTemplate(i)
        hydrogen_indices, heavy_indices = _process_unmapped_atoms(product)
        unmapped_hydrogen_indices[1].append(hydrogen_indices)
        unmapped_heavy_indices[1].append(heavy_indices)

    # ---------------------- POPULATE STATIC HYDROGENS ---------------------
    # Initialize static_hydrogen_indices with the same structure as other indices lists
    static_hydrogen_indices = ([], [])

    # Iterate through reactants and add to the first list in the tuples
    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        static_hydrogens = _process_static_hydrogens(reactant, center_atom_indices[0][i])
        static_hydrogen_indices[0].append(static_hydrogens)

    # Iterate through products and add to the second list in the tuples
    for i in range(reaction.GetNumProductTemplates()):
        product = reaction.GetProductTemplate(i)
        static_hydrogens = _process_static_hydrogens(product, center_atom_indices[1][i])
        static_hydrogen_indices[1].append(static_hydrogens)

    # ==================================================================
    #                       EXTRACT THE OPERATOR
    # ==================================================================

    # Start by duplicating center_atom_indices as keep_atom_indices
    keep_atom_indices = ([], [])

    for i in range(len(center_atom_indices[0])):
        keep_atom_indices[0].append(set(center_atom_indices[0][i]))

    for i in range(len(center_atom_indices[1])):
        keep_atom_indices[1].append(set(center_atom_indices[1][i]))

    # Add or remove indices according to the input booleans specified
    if include_sigma:
        for i in range(len(sigma_atom_indices[0])):
            keep_atom_indices[0][i].update(sigma_atom_indices[0][i])
        for i in range(len(sigma_atom_indices[1])):
            keep_atom_indices[1][i].update(sigma_atom_indices[1][i])

    if include_pi:
        for i in range(len(pi_atom_indices[0])):
            keep_atom_indices[0][i].update(pi_atom_indices[0][i])
        for i in range(len(pi_atom_indices[1])):
            keep_atom_indices[1][i].update(pi_atom_indices[1][i])

    if include_unmapped_hydrogens:
        for i in range(len(unmapped_hydrogen_indices[0])):
            keep_atom_indices[0][i].update(unmapped_hydrogen_indices[0][i])
        for i in range(len(unmapped_hydrogen_indices[1])):
            keep_atom_indices[1][i].update(unmapped_hydrogen_indices[1][i])

    if include_unmapped_heavy_atoms:
        for i in range(len(unmapped_heavy_indices[0])):
            keep_atom_indices[0][i].update(unmapped_heavy_indices[0][i])
        for i in range(len(unmapped_heavy_indices[1])):
            keep_atom_indices[1][i].update(unmapped_heavy_indices[1][i])

    if not include_static_hydrogens:
        for i in range(len(static_hydrogen_indices[0])):
            keep_atom_indices[0][i].difference_update(static_hydrogen_indices[0][i])
        for i in range(len(static_hydrogen_indices[1])):
            keep_atom_indices[1][i].difference_update(static_hydrogen_indices[1][i])

    # Create remove_bond_indices and remove_atom_indices
    remove_bond_indices = ([], [])
    remove_atom_indices = ([], [])

    # Process reactants
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

    # Process products
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

    # Remove bonds and atoms from the reaction
    new_reaction = rdChemReactions.ChemicalReaction()

    # Process reactants
    for i in range(reaction.GetNumReactantTemplates()):
        reactant = reaction.GetReactantTemplate(i)
        editable_reactant = Chem.EditableMol(reactant)
        for bond_idx in remove_bond_indices[0][i]:
            bond = reactant.GetBondWithIdx(bond_idx)
            editable_reactant.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for atom_idx in sorted(remove_atom_indices[0][i], reverse=True):
            editable_reactant.RemoveAtom(atom_idx)
        new_reaction.AddReactantTemplate(editable_reactant.GetMol())

    # Process products
    for i in range(reaction.GetNumProductTemplates()):
        product = reaction.GetProductTemplate(i)
        editable_product = Chem.EditableMol(product)
        for bond_idx in remove_bond_indices[1][i]:
            bond = product.GetBondWithIdx(bond_idx)
            editable_product.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for atom_idx in sorted(remove_atom_indices[1][i], reverse=True):
            editable_product.RemoveAtom(atom_idx)
        new_reaction.AddProductTemplate(editable_product.GetMol())

    # Convert to SMIRKS
    smirks = rdChemReactions.ReactionToSmarts(new_reaction)
    return smirks

