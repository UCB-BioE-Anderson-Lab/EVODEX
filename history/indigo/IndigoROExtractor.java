package org.ucb.c5.roextractor.indigo;

import com.epam.indigo.Indigo;
import com.epam.indigo.IndigoInchi;
import com.epam.indigo.IndigoObject;
import com.epam.indigo.IndigoRenderer;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.ucb.c5.Pair;

/**
 *
 * @author Mindi
 * @author J. Christopher Anderson
 */
public class IndigoROExtractor {

    public static final int ELECTRONIC = 1;
    public static final int CORE = 2;
    public static final int VALENCE = 3;
    public static final int BONDING = 4;

    public static boolean verbose = false;

    public String run(IndigoObject mappedRxn, boolean stereo, boolean hydrogenated, boolean complete, int rotype) {
        IndigoObject mrxn = mappedRxn.clone();

        if (verbose) {
            saveRxnImg(mrxn, "IndigoROExtractor-clone.svg");
        }

        //If stereo is true, then the stereochemistry will be stripped from the molecule
        if (stereo) {
            for (IndigoObject mol : mrxn.iterateMolecules()) {
                mol.clearStereocenters();
            }
        }
        
        if (verbose) {
            saveRxnImg(mrxn, "IndigoROExtractor-postStereo.svg");
        }

        //Index out all the atoms and molecules
        IndigoObject[] molObjects = new IndigoObject[mrxn.countMolecules()];
        IndigoObject[][] atomObjects = new IndigoObject[mrxn.countMolecules()][];
        for (int m = 0; m < mrxn.countMolecules(); m++) {
            IndigoObject mol = mrxn.getMolecule(m);
            molObjects[m] = mol;
            atomObjects[m] = new IndigoObject[mol.countAtoms()];
            for (int a = 0; a < mol.countAtoms(); a++) {
                IndigoObject atom = mol.getAtom(a);
                atomObjects[m][a] = atom;
            }
        }

        //Index out all the bonds
        IndigoObject[][] bondObjects = new IndigoObject[mrxn.countMolecules()][];
        for (int m = 0; m < mrxn.countMolecules(); m++) {
            IndigoObject mol = mrxn.getMolecule(m);
            bondObjects[m] = new IndigoObject[mol.countBonds()];
            for (int b = 0; b < mol.countBonds(); b++) {
                IndigoObject bond = mol.getBond(b);
                bondObjects[m][b] = bond;
            }
        }

        //Create arrays to hold the atom reaction centers
        Set<Pair<Integer, Integer>> centerAtoms = new HashSet<>();
        Set<Pair<Integer, Integer>> centerBonds = new HashSet<>();

        //Create arrays to hold the atom indices to be retained
        Set<Pair<Integer, Integer>> keepAtoms = new HashSet<>();
        Set<Pair<Integer, Integer>> keepBonds = new HashSet<>();

        //Identify any exchangeable hydrogens
        Set<Integer> unmappers = new HashSet<>();
        for (IndigoObject mol : mrxn.iterateMolecules()) {
            for (int a = 0; a < mol.countAtoms(); a++) {
                IndigoObject atom = mol.getAtom(a);

                //Only deal with hydrogens
                if (atom.atomicNumber() != 85) {
                    continue;
                }

                //Gather the hydrogens next to heteroatoms
                for (IndigoObject nei : atom.iterateNeighbors()) {
                    IndigoObject bond = nei.bond();
                    IndigoObject neiAtom = atomObjects[mol.index()][nei.index()];
                    int atno = neiAtom.atomicNumber();
                    if (atno == 7 || atno == 8 || atno == 16) {
                        unmappers.add(mrxn.atomMappingNumber(atom));
                    }
                }
            }
        }

        //Unmap the exchangeable Hydrogens
        for (IndigoObject mol : mrxn.iterateMolecules()) {
            for (int a = 0; a < mol.countAtoms(); a++) {
                IndigoObject atom = mol.getAtom(a);
                if (unmappers.contains(mrxn.atomMappingNumber(atom))) {
                    mrxn.setAtomMappingNumber(atom, 0);
                }
            }
        }

        if (verbose) {
            saveRxnImg(mrxn, "IndigoROExtractor-postExchangeH.svg");
        }

        if (verbose) {
            System.out.println("\n*********\nAll Bonds\n\t\t\tFrom:\t\t\tTo:");
            System.out.println("m\tb\tcode\tindex\tatno\t\tindex\tatno");
        }

        //Iterate through the bonds to identify rxn centers
        for (int m = 0; m < bondObjects.length; m++) {
            for (int b = 0; b < bondObjects[m].length; b++) {
                IndigoObject bond = bondObjects[m][b];

                //Get the enum that describes the mapping of each reaction center
                int code = mrxn.reactingCenter(bond);

                if (verbose) {
                    System.out.print(m + "\t" + b + "\t" + code + "\t");
                    System.out.print(bond.source().index() + "\t" + bond.source().atomicNumber());
                    System.out.print("\t\t");
                    System.out.println(bond.destination().index() + "\t" + bond.destination().atomicNumber());
                }

                if (code == Indigo.RC_UNCHANGED) {    //Atoms on both sides, no bonds changed
                    continue;
                }
                if (code == Indigo.RC_CENTER) {       //Not sure what that means, but isn't encountered much if ever
                    continue;
                }

                IndigoObject atomOne = bond.source();
                IndigoObject atomTwo = bond.destination();

                if (complete) {
                    //With the complete option, all changed atoms get included
                    keepBonds.add(new Pair(m, b));
                    keepAtoms.add(new Pair(m, atomOne.index()));
                    keepAtoms.add(new Pair(m, atomTwo.index()));
                    centerAtoms.add(new Pair(m, atomOne.index()));
                    centerAtoms.add(new Pair(m, atomTwo.index()));
                    centerBonds.add(new Pair(m, b));
                } else {
                    //Matched scenario, only include the atoms and bonds that are not atommap == 0
                    boolean onekeep = false;
                    boolean twokeep = false;
                    if (mrxn.atomMappingNumber(atomOne) != 0) {
                        keepAtoms.add(new Pair(m, atomOne.index()));
                        centerAtoms.add(new Pair(m, atomOne.index()));
                        onekeep = true;
                    }
                    if (mrxn.atomMappingNumber(atomTwo) != 0) {
                        keepAtoms.add(new Pair(m, atomTwo.index()));
                        centerAtoms.add(new Pair(m, atomTwo.index()));
                        twokeep = true;
                    }
                    if (onekeep && twokeep) {
                        keepBonds.add(new Pair(m, bond.index()));
                        centerBonds.add(new Pair(m, b));
                    }
                }
            }
        }

        if (verbose) {
            /**
             * m is the molecule, substrate = 0, product = 1 a is the atom index
             * within that molecule atno is the atomic number (85 is At, C=6,
             * O=8) atmap is the assigned atom mapping number
             */
            System.out.println("\n*********\nKeep Atoms after collecting reaction centers");
            System.out.println("m\ta\tatno\tatmap");

            for (Pair<Integer, Integer> apair : keepAtoms) {
                IndigoObject atom = atomObjects[apair.getKey()][apair.getValue()];
                System.out.println(apair.getKey() + "\t" + apair.getValue() + "\t" + atom.atomicNumber() + "\t" + mrxn.atomMappingNumber(atom));
            }
        }

        if (rotype == ELECTRONIC) {
            //For each reaction center Atom, add the sigma-attached bonds and atoms
            for (Pair<Integer, Integer> index : centerAtoms) {
                IndigoObject croAtom = atomObjects[index.getKey()][index.getValue()];
                for (IndigoObject nei : croAtom.iterateNeighbors()) {
                    IndigoObject bond = nei.bond();
                    keepBonds.add(new Pair(index.getKey(), bond.index()));
                    keepAtoms.add(new Pair(index.getKey(), bond.source().index()));
                    keepAtoms.add(new Pair(index.getKey(), bond.destination().index()));
                }
            }

            //For each keepAtom, add anything in conjugation
            Set<Pair<Integer, Integer>> workList = new HashSet<>();
            workList.addAll(keepAtoms);

            for (IndigoObject mol : mrxn.iterateMolecules()) {
                for (int a = 0; a < mol.countAtoms(); a++) {
                    IndigoObject atom = mol.getAtom(a);
                    if (!keepAtoms.contains(new Pair(mol.index(), a))) {
                        continue;
                    }
                    IndigoObject keeper = mol.getAtom(a);
                    if (isSP3(keeper)) {
                        continue;
                    }
                    keepBonds.addAll(addConjugated(mol, keeper, keepAtoms));
                }
            }
        }

        if (verbose) {
            System.out.println("\n*********\nAfter Electronic block");
            System.out.println("m\ta\tatno\tatmap");

            for (Pair<Integer, Integer> apair : keepAtoms) {
                IndigoObject atom = atomObjects[apair.getKey()][apair.getValue()];
                System.out.println(apair.getKey() + "\t" + apair.getValue() + "\t" + atom.atomicNumber() + "\t" + mrxn.atomMappingNumber(atom));
            }
        }

        if (hydrogenated) {
            //Remove certain hydrogens that are mapped
            Set<Pair<Integer, Integer>> keepcopy = new HashSet<>(keepAtoms);
            for (Pair<Integer, Integer> keeper : keepcopy) {
                int m = keeper.getKey();
                int a = keeper.getValue();
                IndigoObject atom = atomObjects[m][a];

                //Only deal with hydrogens
                int atno = atom.atomicNumber();
                if (atno != 85) {
                    continue;
                }

                //Exclude unmapped
                int atmap = mrxn.atomMappingNumber(atom);
                if (atmap == 0) {
                    continue;
                }

                //Exclude reactive center hydrogens
                if (centerAtoms.contains(new Pair(m, a))) {
                    continue;
                }

                keepAtoms.remove(keeper);
                for (IndigoObject nei : atom.iterateNeighbors()) {
                    IndigoObject bond = nei.bond();
                    keepBonds.remove(new Pair(m, bond.index()));
                }
            }
        }

        if (verbose) {
            System.out.println("\n*********\nAfter Hydrogenize block");
            System.out.println("m\ta\tatno\tatmap");

            for (Pair<Integer, Integer> apair : keepAtoms) {
                IndigoObject atom = atomObjects[apair.getKey()][apair.getValue()];
                System.out.println(apair.getKey() + "\t" + apair.getValue() + "\t" + atom.atomicNumber() + "\t" + mrxn.atomMappingNumber(atom));
            }
        }

        //Gather up bonds that should be tossed
        Set<IndigoObject> tossBonds = new HashSet<>();
        for (int m = 0; m < mrxn.countMolecules(); m++) {
            IndigoObject mol = mrxn.getMolecule(m);
            for (int b = 0; b < mol.countBonds(); b++) {
                IndigoObject bond = bondObjects[m][b];
                tossBonds.add(bond);
            }
        }
        for (Pair<Integer, Integer> indices : keepBonds) {
            IndigoObject keepbond = bondObjects[indices.getKey()][indices.getValue()];
            tossBonds.remove(keepbond);
        }

        //Gather up atoms that should be tossed
        Set<IndigoObject> tossAtoms = new HashSet<>();
        for (int m = 0; m < mrxn.countMolecules(); m++) {
            IndigoObject mol = mrxn.getMolecule(m);
            for (int a = 0; a < mol.countAtoms(); a++) {
                IndigoObject atom = atomObjects[m][a];
                tossAtoms.add(atom);
            }
        }
        for (Pair<Integer, Integer> indices : keepAtoms) {
            IndigoObject keepatom = atomObjects[indices.getKey()][indices.getValue()];
            tossAtoms.remove(keepatom);
        }

        //Convert Astatines back to Hydrogens
        for (IndigoObject mol : mrxn.iterateMolecules()) {
            for (int i = 0; i < mol.countAtoms(); i++) {
                IndigoObject atom = mol.getAtom(i);
                if (atom.atomicNumber() == 85) {
                    atom.resetAtom("H");
                }
            }
        }

        //Remove the tossed bonds
        for (IndigoObject bond : tossBonds) {
            bond.remove();
        }

        //Remove the tossed atoms
        for (IndigoObject atom : tossAtoms) {
            atom.remove();
        }
        
        if(verbose) {
            saveRxnImg(mrxn, "IndigoROExtractor_postToss.svg");
        }

        /**
         * If an RO is of the valence type, up to here it has behaved as core.
         * None of the electronic block would be called. It could still be
         * matched or complete, and it could be hydrogenized or not.
         *
         * Regardless, the VALENCE RO is obtained by eliminating the identity of
         * all the atoms. Here that is done by converting to a *.
         */
        if (rotype == VALENCE) {
            for (IndigoObject mol : mrxn.iterateMolecules()) {
                for (IndigoObject atom : mol.iterateAtoms()) {
                    atom.resetAtom("*");
                }
            }
        }

        //Remove the extra stuff at the end of the smirks and return
        String smirks = mrxn.smiles();
        smirks = smirks.split("\\|")[0];
        return smirks.trim();
    }

    /**
     * Iterates through the bonds on the input atom. If any of those bonds are
     * double or triple bonds, it returns false. Otherwise, it is SP3.
     *
     * @param atom
     * @return
     */
    private boolean isSP3(IndigoObject atom) {
        for (IndigoObject nei : atom.iterateNeighbors()) {
            IndigoObject bond = nei.bond();
            int order = bond.bondOrder();
            if (order > 1) {
                return false;
            }
        }
        return true;

    }

    private Set<Pair<Integer, Integer>> addConjugated(IndigoObject mol, IndigoObject croAtom, Set<Pair<Integer, Integer>> keepAtoms) {
        Set<Pair<Integer, Integer>> out = new HashSet<>();

        for (IndigoObject nei : croAtom.iterateNeighbors()) {
            int neiIndex = nei.index();   // index of neighbor atom
            IndigoObject bond = nei.bond();
            IndigoObject neiAtom = mol.getAtom(neiIndex);
            if (!isSP3(neiAtom)) {
                out.add(new Pair(mol.index(), bond.index()));
                if (!keepAtoms.contains(new Pair(mol.index(), neiAtom.index()))) {
                    keepAtoms.add(new Pair(mol.index(), neiAtom.index()));
                    out.addAll(addConjugated(mol, neiAtom, keepAtoms));
                }
            }
        }
        return out;
    }

    private void saveRxnImg(IndigoObject mrxn, String path) {
        Indigo indigo = new Indigo();
        IndigoRenderer renderer = new IndigoRenderer(indigo);
        indigo.setOption("render-output-format", "svg");
        indigo.setOption("render-implicit-hydrogens-visible", false);
        indigo.setOption("embedding-uniqueness", "none");
        String smiles = mrxn.smiles();
        IndigoObject reimport = indigo.loadReaction(smiles);
        renderer.renderToFile(reimport, path);
    }

    public static void main(String[] args) throws Exception {
        verbose = true;

        // NOT_CENTER:
        System.out.println("RC_NOT_CENTER " + Indigo.RC_NOT_CENTER);
        // UNMARKED: no atom map because atom is only on one side of the reaction
        System.out.println("RC_UNMARKED " + Indigo.RC_UNMARKED);
        // CENTER:
        System.out.println("RC_CENTER " + Indigo.RC_CENTER);
        // UNCHANGED: the bond is the same on both sides of the equation
        System.out.println("RC_UNCHANGED " + Indigo.RC_UNCHANGED);
        // MADE_OR_BROKEN: the bond is either made or broken during the reaction
        System.out.println("RC_MADE_OR_BROKEN " + Indigo.RC_MADE_OR_BROKEN);
        // ORDER_CHANGED: went from double bond to single bond, etc
        System.out.println("RC_ORDER_CHANGED " + Indigo.RC_ORDER_CHANGED);

//        String reaction = "[C:11]([H:12])([H:13])([H:14])[C:7]([H:8])=[C:6]([H:10])[C:1]([H:3])([H])[O:5]([H])>>[C:11]([H:12])([H:13])([H:14])[C:7]([H:8])=[C:6]([H:10])[C:1]([H:3])=[O:5]"; // proponyl alcohol to methyl acrolein
//        System.out.println(newReaction);
//        String reaction = "[C:7]([H:8])([H:9])=[C:6]([H:10])[C:1]([H:3])([H])[O:5]([H])>>[C:7]([H:8])([H:9])=[C:6]([H:10])[C:1]([H:3])=[O:5]"; // allyl alcohol to acrolein
//        String reaction = "[N:1]([H:2])[H] >> [N:1]([H:2])C(=O)C"; // acetylation
//        String newReaction = reaction.replaceAll("\\[H", "[At");
//        String reaction = "[C:3](/[C:1]([O:20][At])([At:14])[At:11])(=[C:4](\\[c:2]1[c:12]([At:16])[c:10]([At:17])[c:6]([O:5][At:8])[c:13]([At:19])[c:15]1[At:18])/[At:7])\\[At:9]>>C(C([O:20][C:1](/[C:3](=[C:4](/[c:2]1[c:12]([At:16])[c:10]([At:17])[c:6]([O:5][At:8])[c:13]([At:19])[c:15]1[At:18])\\[At:7])/[At:9])([At:14])[At:11])=O)([At])([At])[At]";
//        String reaction = "[c:1]1([c:2]([C:12](=[O:14])C(O[At:13])=O)[c:5]([At:7])[c:3]([At:8])[c:4]([O:6][At:9])[c:1]1[At:10])[At:11]>>[c:1]1([c:2]([C:12](=[O:14])[At:13])[c:5]([At:7])[c:3]([At:8])[c:4]([O:6][At:9])[c:1]1[At:11])[At:10]";
        String reaction = "[C:1]([C:1]([C:2]([C:3]([C:4]([C:5]([C:6]([C:7]([C:8]([C:9]([C:10]([C:11]([C:12]([C:13]([C:14]([C:46](=[O:47])[At])([At:15])[At:29])([At:17])[At:31])([At:30])[At:16])([At:33])[At:19])([At:18])[At:32])([At:21])[At:35])([At:20])[At:34])([At:37])[At:23])([At:22])[At:36])([At:39])[At:25])([At:24])[At:38])([At:41])[At:27])([At:26])[At:40])([At:28])[At:42])([At:43])([At:45])[At:44]>>[C:1]([C:1]([C:2]([C:3]([C:4]([C:5]([C:6]([C:7]([C:8]([C:9]([C:10]([C:11]([C:12]([C:13]([C:14]([C:46](SC(C(N=C(O[At])C(C(N=C(O[At])C(O[At])(C(C(OP(OP(OC([C@]1(O[C@]([n]2c3[n]c([n]c(N([At])[At])c3[n]c2[At])[At])([At])[C@@](O[At])([At])[C@@]1(OP(=O)(O[At])O[At])[At])[At])([At])[At])(=O)O[At])(=O)O[At])([At])[At])(C([At])([At])[At])C([At])([At])[At])[At])([At])[At])([At])[At])([At])[At])([At])[At])=[O:47])([At:29])[At:15])([At:31])[At:17])([At:30])[At:16])([At:33])[At:19])([At:18])[At:32])([At:21])[At:35])([At:20])[At:34])([At:23])[At:37])([At:36])[At:22])([At:39])[At:25])([At:38])[At:24])([At:41])[At:27])([At:40])[At:26])([At:44])[At:42])([At:28])([At:43])[At:45] |w:78,a:89,91,107,111|";
        // creating the rxn molecule
        Indigo indigo = new Indigo();

        IndigoObject mappedRxn = indigo.loadReaction(reaction);
        mappedRxn.aromatize();
        mappedRxn.correctReactingCenters();

        IndigoRenderer renderer = new IndigoRenderer(indigo);
        renderer.renderToFile(mappedRxn, "IndigoROExtractor-mapped.svg");

        IndigoObject aclone = mappedRxn.clone();
        for (IndigoObject mol : aclone.iterateMolecules()) {
            for (int i = 0; i < mol.countAtoms(); i++) {
                IndigoObject atom = mol.getAtom(i);
                aclone.setAtomMappingNumber(atom, i);
            }
        }
        renderer.renderToFile(aclone, "IndigoROExtractor-indices.svg");

        IndigoROExtractor roe = new IndigoROExtractor();

        // mappedRxn,  stereo,  hydrogenated,  complete,  rotype
        String ro = roe.run(mappedRxn, false, true, true, IndigoROExtractor.ELECTRONIC);

        System.out.println("\n*********\nResult\n" + ro);
        IndigoObject rxn1 = indigo.loadQueryReaction(ro);
        renderer.renderToFile(rxn1, "IndigoROExtractor-finalresult.svg");
    }

}
