package org.ucb.c5.roextractor.chemaxon;

import chemaxon.formats.MolImporter;
import chemaxon.standardizer.Standardizer;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.RxnMolecule;
import java.io.File;
import java.util.HashSet;
import java.util.Set;
import org.ucb.c5.ro.utils.ChemAxonUtils;

/**
* This Function inputs a concrete reaction (the variable mapped) that has been
* cleaned by standardization. If it is a "half" reaction, then the cofactors
* should be removed prior to invoking this. Also, the atoms should already have
* been mapped by SkeletonMapper, ChangeMapper, or otherwise. It returns a reaction 
* operator of the chosen configuration as a RxnMolecule object.
*
* @author J. Christopher Anderson
*/
public class ROExtractor {
   
   public static final int ELECTRONIC = 1;
   public static final int CORE = 2;
   public static final int VALENCE = 3;
   public static final int BONDING = 4;

   /*
"cofactors"  Half/Full   (include cofactors or not)
"complete" Complete/Matched (include all atoms, or just the reacting set)
"hydrogenated" Hydrogenate/Not  (include additional hydrogens to handle valency)
Electronic/Core/Bond/Valency  (include just the core
"stereo"  stereochemistry (include, or remove)
   
   I think this comes out to 24 possible ROs
    */
   /**
    * @param mapped The standardized and mapped, concrete RxnMolecule to be
    * extracted
    * @param stereo Retain the stereochemistry or make the RO achiral
    * @param complete Retain unmapped atoms or include only matching atoms
    * @param hydrogenated Retain hydrogens to preserve valency, or omit them
    * @param rotype The type of reaction operator (electronic:1, core:2, valence:3, or
    * bonding:4)
    * @return the reaction operator of the specified configuration
    */
   public RxnMolecule run(RxnMolecule mappedRxn, boolean stereo, boolean hydrogenated, boolean complete, int rotype) {
       //Make a copy as the returned object
       RxnMolecule mapped = mappedRxn.clone();

       //Put the hybridization state on each atom
       mapped.calcHybridization();

       //Gather up the map indices of atoms that are reaction centers
       Set<MolBond> keepBonds = new HashSet<>();
       Set<Integer> rxnCenters = new HashSet<>();
       for (MolBond bond : mapped.getBondArray()) {
           if ((bond.getFlags() & MolBond.REACTING_CENTER_MASK) == 0) {
               continue;
           }
           keepBonds.add(bond);
           MolAtom one = bond.getAtom1();
           rxnCenters.add(one.getAtomMap());
           MolAtom two = bond.getAtom2();
           rxnCenters.add(two.getAtomMap());
       }

       //Gather up all the atoms that aren't reaction centers
       Set<MolAtom> keepAtoms = new HashSet<>();
       Set<MolAtom> changers = new HashSet<>();
       for (int i = 0; i < mapped.getAtomCount(); i++) {
           MolAtom atom = mapped.getAtom(i);
           if (rxnCenters.contains(atom.getAtomMap())) {
               keepAtoms.add(atom);
               changers.add(atom);
           }
       }

       if (rotype == ELECTRONIC) {
           //For each changer, add the sigma-attached bonds and atoms
           for (MolAtom croAtom : changers) {
               for (int i = 0; i < croAtom.getBondCount(); i++) {
                   MolBond bond = croAtom.getBond(i);
                   keepBonds.add(bond);
                   keepAtoms.add(bond.getOtherAtom(croAtom));
               }
           }

           //For each keepAtom, add anything in conjugation
           Set<MolAtom> workList = new HashSet<>();
           workList.addAll(keepAtoms);

           for (MolAtom keeper : workList) {
               int sp = keeper.getHybridizationState();
               if (sp == MolAtom.HS_SP2 || sp == MolAtom.HS_SP) {
                   keepBonds.addAll(addConjugated(keeper, keepAtoms));
               }
           }
       }

       //Gather up bonds that should be tossed
       Set<MolBond> tossBonds = new HashSet<>();
       for (int i = 0; i < mapped.getBondCount(); i++) {
           MolBond bond = mapped.getBond(i);
           if (!keepBonds.contains(bond)) {
               tossBonds.add(bond);
           }
       }

       //Gather up atoms that should be tossed
       Set<MolAtom> tossAtoms = new HashSet<>();
       for (int i = 0; i < mapped.getAtomCount(); i++) {
           MolAtom atom = mapped.getAtom(i);
           if (!keepAtoms.contains(atom)) {
               tossAtoms.add(atom);
           }
       }

       if (!hydrogenated) {
           //Remove any hydrogens that are mapped
           for (MolAtom atom : mapped.getAtomArray()) {
               //Only deal with hydrogens
               if (atom.getAtno() > 1) {
                   continue;
               }
               //Exclude unmapped
               if (atom.getAtomMap() == 0) {
                   continue;
               }
               tossAtoms.add(atom);
               MolBond bond = atom.getBond(0);
               tossBonds.add(bond);
           }
       }

       //Remove the tossed bonds
       for (MolBond bond : tossBonds) {
           mapped.removeBond(bond);
       }

       //Remove the tossed atoms
       for (MolAtom tossme : tossAtoms) {
           mapped.removeAtom(tossme);
       }

       if (rotype == VALENCE) {
           //Make everything single bonds
           for (MolBond bond : mapped.getBondArray()) {
               bond.setType(1);
           }
       }

       if (rotype == BONDING) {
           //Remove all the atoms
           Set<MolAtom> finalToss = new HashSet<>();
           for (MolAtom atom : mapped.getAtomArray()) {
               finalToss.add(atom);
           }
           for (MolAtom atom : finalToss) {
               mapped.removeAtom(atom);
           }
       }

       if (!stereo) {
           //Remove any chirality
           Standardizer std1 = new Standardizer("clearstereo");
           std1.standardize(mapped);
       }

       return mapped;
   }

   private Set<MolBond> addConjugated(MolAtom croAtom, Set<MolAtom> keepers) {
       Set<MolBond> out = new HashSet<>();

       for (int i = 0; i < croAtom.getBondCount(); i++) {
           MolBond bond = croAtom.getBond(i);
           MolAtom nei = bond.getOtherAtom(croAtom);

           int sp = nei.getHybridizationState();
           if (sp == MolAtom.HS_SP2 || sp == MolAtom.HS_SP) {
               out.add(bond);
               if (!keepers.contains(nei)) {
                   keepers.add(nei);
                   out.addAll(addConjugated(nei, keepers));
               }
           }
       }
       return out;
   }

   public static void main(String[] args) throws Exception {
       ChemAxonUtils.license();

       File output = new File("output");
       if (!output.exists()) {
           output.mkdir();
       }

//        String reaction = "CCCCO>>CCCC=O"; //alcohol to aldehyde
//        String reaction = "OCC(OP(O)(O)=O)C(O)=O>>OC(COP(O)(O)=O)C(O)=O"; //2-PG >> 3-PG
//        String reaction = "CCCC=CC(=O)N>>CCCC=CC(=O)O"; //amide to acid
//        String reaction = "CCCc1ccccc1C(=O)N>>CCCc1ccccc1C(=O)O"; //amide to acid
//        String reaction = "C[C@@H](C(=O)O)O>>CC(=O)C(=O)O"; //lactate to pyruvate (lactate dehydrogenase)

       //acrylate + coenzyme A + H+ → acryloyl-CoA + H2O
       String reaction = "CCCC=CC(=O)O"; 
       reaction += ">>";
       reaction += "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)[C@H](C(=O)NCCC(=O)NCCSC(=O)C=CCCC)O";
       
       RxnMolecule rxn = RxnMolecule.getReaction(MolImporter.importMol(reaction));

       Standardizer std = new Standardizer("addexplicitH");
       std.standardize(rxn);

       RxnMolecule mapped = new SkeletonMapper().run(rxn);
       ChemAxonUtils.savePNGImage(mapped, "output/example_mapped.png");

       ROExtractor extractor = new ROExtractor();

       RxnMolecule shchERO = extractor.run(mapped, true, true, true, ROExtractor.ELECTRONIC);
       ChemAxonUtils.saveSMIRKSPNGImage(shchERO, "output/example_shchERO.png");

       RxnMolecule hchERO = extractor.run(mapped, false, true, true, ROExtractor.ELECTRONIC);
       ChemAxonUtils.saveSMIRKSPNGImage(hchERO, "output/example_hchERO.png");

       RxnMolecule schERO = extractor.run(mapped, true, true, false, ROExtractor.ELECTRONIC);
       ChemAxonUtils.saveSMIRKSPNGImage(schERO, "output/example_schERO.png");

       RxnMolecule shchCRO = extractor.run(mapped, true, true, true, ROExtractor.CORE);
       ChemAxonUtils.saveSMIRKSPNGImage(shchCRO, "output/example_shchCRO.png");

//        RxnMolecule hmhCRO = extractor.run(mapped, false, false, true, ROExtractor.CORE);
//        ChemAxonUtils.saveSMIRKSPNGImage(hmhCRO, "output/example_hmhCRO.png");
   }
}
