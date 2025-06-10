import pytest
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from evodex.decofactor import remove_cofactors

@pytest.fixture(scope="module")
def load_reactions_with_astatine():
    reactions = []
    with open('tests/data/splitting_test_data.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            reactions.append(row['atom_mapped'])
    return reactions

def test_remove_cofactors(load_reactions_with_astatine):
    for reaction_smiles in load_reactions_with_astatine:
        try:
            partial_reaction = remove_cofactors(reaction_smiles)
            
            # Ensure the partial reaction is not None or empty
            assert partial_reaction is not None, "No partial reactions generated"
            
            # Check if the partial reaction is ">>"
            if partial_reaction == ">>":
                continue  # Accept ">>" as valid behavior
            
            try:
                rxn = rdChemReactions.ReactionFromSmarts(partial_reaction, useSmiles=True)
                
                # Check that each partial reaction can be parsed as a reaction SMILES
                assert rxn is not None, f"Invalid reaction SMILES: {partial_reaction}"
                
                # Ensure all molecules contain at least one mapped atom on both sides
                substrate_atom_maps = set()
                for mol in rxn.GetReactants():
                    for atom in mol.GetAtoms():
                        if atom.GetAtomMapNum() > 0:
                            substrate_atom_maps.add(atom.GetAtomMapNum())
                
                product_atom_maps = set()
                for mol in rxn.GetProducts():
                    for atom in mol.GetAtoms():
                        if atom.GetAtomMapNum() > 0:
                            product_atom_maps.add(atom.GetAtomMapNum())
                
                assert len(substrate_atom_maps) > 0, f"No mapped atoms in reactants: {partial_reaction}"
                assert len(product_atom_maps) > 0, f"No mapped atoms in products: {partial_reaction}"
                
                # Ensure all atom maps are unique on either side
                assert len(substrate_atom_maps) == sum(1 for mol in rxn.GetReactants() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0), f"Duplicate atom maps in reactants: {partial_reaction}"
                assert len(product_atom_maps) == sum(1 for mol in rxn.GetProducts() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0), f"Duplicate atom maps in products: {partial_reaction}"
                
                # Ensure the set of maps present in the substrates matches those in the products
                assert substrate_atom_maps == product_atom_maps, f"Mismatch between reactant and product atom maps: {partial_reaction}"
            
            except Exception as e:
                pytest.fail(f"Failed to validate partial reaction: {partial_reaction}, Error: {e}")
        
        except Exception as e:
            pytest.fail(f"remove_cofactors raised an error: {e}")

if __name__ == "__main__":
    pytest.main()


# Test for contains_cofactor
def test_contains_cofactor_classification():
    from evodex.decofactor import contains_cofactor

    # Placeholder reactions containing cofactors (edit these with real examples)
    reactions_with_cofactors = [
        # Full oxidoreductase
        "[C:1]([C:2]([C:3]([C:4]([C:5]([C:6]([C:7]([C:8]([C:9]([C:10]([C:11]([C:12]([C:13]([C:14]([C:15]([C:16]([O:17][C:18]([C:19](=[O:20])[C:21]([O:22][P:23](=[O:24])([O:25][H:71])[O:26][H:72])([H:73])[H:74])([H:75])[H:76])([H:77])[H:78])([H:79])[H:80])([H:81])[H:82])([H:83])[H:84])([H:85])[H:86])([H:87])[H:88])([H:89])[H:90])([H:91])[H:92])([H:93])[H:94])([H:95])[H:96])([H:97])[H:98])([H:99])[H:100])([H:101])[H:102])([H:103])[H:104])([H:105])[H:106])([H:107])([H:108])[H:109].[HH].[H][C:70]1([H:135])[C:30]([C:28]([N:27]([H:137])[H:138])=[O:29])=[C:31]([H:110])[N:32]([C@:33]2([H:111])[O:34][C@@:35]([C:36]([O:37][P:38](=[O:39])([O:40][H:112])[O:41][P:42](=[O:43])([O:44][H:113])[O:45][C:46]([C@@:47]3([H:114])[O:48][C@:49]([N:50]4:[C:51]([H:115]):[N:52]:[C:53]5:[C:54]([N:55]([H:116])[H:117]):[N:56]:[C:57]([H:118]):[N:58]:[C:59]:4:5)([H:119])[C@@:60]([O:61][H:120])([H:121])[C@@:62]3([O:63][H:122])[H:123])([H:124])[H:125])([H:126])[H:127])([H:128])[C@:64]([O:65][H:129])([H:130])[C@:66]2([O:67][H:131])[H:132])[C:68]([H:133])=[C:69]1[H:134]>>[H][C@:19]([C:18]([O:17][C:16]([C:15]([C:14]([C:13]([C:12]([C:11]([C:10]([C:9]([C:8]([C:7]([C:6]([C:5]([C:4]([C:3]([C:2]([C:1]([H:107])([H:108])[H:109])([H:105])[H:106])([H:103])[H:104])([H:101])[H:102])([H:99])[H:100])([H:97])[H:98])([H:95])[H:96])([H:93])[H:94])([H:91])[H:92])([H:89])[H:90])([H:87])[H:88])([H:85])[H:86])([H:83])[H:84])([H:81])[H:82])([H:79])[H:80])([H:77])[H:78])([H:75])[H:76])([O:20][H])[C:21]([O:22][P:23](=[O:24])([O:25][H:71])[O:26][H:72])([H:73])[H:74].[N:27]([C:28](=[O:29])[C:30]1:[C:31]([H:110]):[N+:32]([C@:33]2([H:111])[O:34][C@@:35]([C:36]([O:37][P:38](=[O:39])([O:40][H:112])[O:41][P:42](=[O:43])([O:44][H:113])[O:45][C:46]([C@@:47]3([H:114])[O:48][C@:49]([N:50]4:[C:51]([H:115]):[N:52]:[C:53]5:[C:54]([N:55]([H:116])[H:117]):[N:56]:[C:57]([H:118]):[N:58]:[C:59]:4:5)([H:119])[C@@:60]([O:61][H:120])([H:121])[C@@:62]3([O:63][H:122])[H:123])([H:124])[H:125])([H:126])[H:127])([H:128])[C@:64]([O:65][H:129])([H:130])[C@:66]2([O:67][H:131])[H:132]):[C:68]([H:133]):[C:69]([H:134]):[C:70]:1[H:135])([H:137])[H:138]",
        # Full oxidoreducatase with At's
        "[#6:1](-[#6:2](-[#6:3](-[#6:4](-[#6:5](-[#6:6](-[#6:7](-[#6:8](-[#6:9](-[#6:10](-[#6:11](-[#6:12](-[#6:13](-[#6:14](-[#6:15](-[#6:16](-[#8:17]-[#6:18](-[#6:19](=[#8:20])-[#6:21](-[#8:22]-[#15:23](=[#8:24])(-[#8:25]-[At:71])-[#8:26]-[At:72])(-[At:73])-[At:74])(-[At:75])-[At:76])(-[At:77])-[At:78])(-[At:79])-[At:80])(-[At:81])-[At:82])(-[At:83])-[At:84])(-[At:85])-[At:86])(-[At:87])-[At:88])(-[At:89])-[At:90])(-[At:91])-[At:92])(-[At:93])-[At:94])(-[At:95])-[At:96])(-[At:97])-[At:98])(-[At:99])-[At:100])(-[At:101])-[At:102])(-[At:103])-[At:104])(-[At:105])-[At:106])(-[At:107])(-[At:108])-[At:109].[AtH].[#7:27](-[#6:28](=[#8:29])-[#6:30]1=[#6:31](-[At:110])-[#7:32](-[#6@:33]2(-[At:111])-[#8:34]-[#6@@:35](-[#6:36](-[#8:37]-[#15:38](=[#8:39])(-[#8:40]-[At:112])-[#8:41]-[#15:42](=[#8:43])(-[#8:44]-[At:113])-[#8:45]-[#6:46](-[#6@@:47]3(-[At:114])-[#8:48]-[#6@:49](-[#7:50]4:[#6:51](-[At:115]):[#7:52]:[#6:53]5:[#6:54](-[#7:55](-[At:116])-[At:117]):[#7:56]:[#6:57](-[At:118]):[#7:58]:[#6:59]:4:5)(-[At:119])-[#6@@:60](-[#8:61]-[At:120])(-[At:121])-[#6@@:62]-3(-[#8:63]-[At:122])-[At:123])(-[At:124])-[At:125])(-[At:126])-[At:127])(-[At:128])-[#6@:64](-[#8:65]-[At:129])(-[At:130])-[#6@:66]-2(-[#8:67]-[At:131])-[At:132])-[#6:68](-[At:133])=[#6:69](-[At:134])-[#6:70]-1(-[At:135])-[At])(-[At:137])-[At:138]>>[#6:1](-[#6:2](-[#6:3](-[#6:4](-[#6:5](-[#6:6](-[#6:7](-[#6:8](-[#6:9](-[#6:10](-[#6:11](-[#6:12](-[#6:13](-[#6:14](-[#6:15](-[#6:16](-[#8:17]-[#6:18](-[#6@@:19](-[#8:20]-[At])(-[#6:21](-[#8:22]-[#15:23](=[#8:24])(-[#8:25]-[At:71])-[#8:26]-[At:72])(-[At:73])-[At:74])-[At])(-[At:75])-[At:76])(-[At:77])-[At:78])(-[At:79])-[At:80])(-[At:81])-[At:82])(-[At:83])-[At:84])(-[At:85])-[At:86])(-[At:87])-[At:88])(-[At:89])-[At:90])(-[At:91])-[At:92])(-[At:93])-[At:94])(-[At:95])-[At:96])(-[At:97])-[At:98])(-[At:99])-[At:100])(-[At:101])-[At:102])(-[At:103])-[At:104])(-[At:105])-[At:106])(-[At:107])(-[At:108])-[At:109].[#7:27](-[#6:28](=[#8:29])-[#6:30]1:[#6:31](-[At:110]):[#7+:32](-[#6@:33]2(-[At:111])-[#8:34]-[#6@@:35](-[#6:36](-[#8:37]-[#15:38](=[#8:39])(-[#8:40]-[At:112])-[#8:41]-[#15:42](=[#8:43])(-[#8:44]-[At:113])-[#8:45]-[#6:46](-[#6@@:47]3(-[At:114])-[#8:48]-[#6@:49](-[#7:50]4:[#6:51](-[At:115]):[#7:52]:[#6:53]5:[#6:54](-[#7:55](-[At:116])-[At:117]):[#7:56]:[#6:57](-[At:118]):[#7:58]:[#6:59]:4:5)(-[At:119])-[#6@@:60](-[#8:61]-[At:120])(-[At:121])-[#6@@:62]-3(-[#8:63]-[At:122])-[At:123])(-[At:124])-[At:125])(-[At:126])-[At:127])(-[At:128])-[#6@:64](-[#8:65]-[At:129])(-[At:130])-[#6@:66]-2(-[#8:67]-[At:131])-[At:132]):[#6:68](-[At:133]):[#6:69](-[At:134]):[#6:70]:1-[At:135])(-[At:137])-[At:138]",
        # Phosphate liberation with At's
        "[#7:1](=[#6:2]1:[#7:3](-[At:42]):[#6:4](=[#8:5]):[#6:6]2:[#7:7]:[#6:8](-[At:43]):[#7:9](-[#6@:10]3(-[At:44])-[#8:11]-[#6@@:12](-[#6:13](-[#8:14]-[#15:15](=[#8:16])(-[#8:17]-[At:45])-[#8:18]-[#15:19](=[#8:20])(-[#8:21]-[At:46])-[#8:22]-[#15:23](=[#8:24])(-[#8:25]-[At:47])-[#8:26]-[At:48])(-[At:49])-[At:50])(-[At:51])-[#6@:27](-[#8:28]-[#15:29](=[#8:30])(-[#8:31]-[At:52])-[#8:32]-[#15:33](=[#8:34])(-[#8:35]-[At:53])-[#8:36]-[At:54])(-[At:55])-[#6@:37]-3(-[#8:38]-[At:56])-[At:57]):[#6:39]:2:[#7:40]:1-[At:58])-[At:59]>>[#7:1](=[#6:2]1:[#7:3](-[At:42]):[#6:4](=[#8:5]):[#6:6]2:[#7:7]:[#6:8](-[At:43]):[#7:9](-[#6@:10]3(-[At:44])-[#8:11]-[#6@@:12](-[#6:13](-[#8:14]-[#15:15](=[#8:16])(-[#8:17]-[At:45])-[#8:18]-[#15:19](=[#8:20])(-[#8:21]-[At:46])-[#8]-[At])(-[At:49])-[At:50])(-[At:51])-[#6@:27](-[#8:28]-[#15:29](=[#8:30])(-[#8:31]-[At:52])-[#8:32]-[#15:33](=[#8:34])(-[#8:35]-[At:53])-[#8:36]-[At:54])(-[At:55])-[#6@:37]-3(-[#8:38]-[At:56])-[At:57]):[#6:39]:2:[#7:40]:1-[At:58])-[At:59].[#8:22](-[#15:23](=[#8:24])(-[#8:25]-[At:47])-[#8:26]-[At:48])-[At]",
        # With implicit H2O
        "CN>>C=O.O",
    ]

    # Placeholder reactions without cofactors (edit these with real examples)
    reactions_without_cofactors = [
        # With H's
        "[C:1]([C:2]([C:3]([C:4]([H:54])([H:55])[H:56])=[O:5])([H:57])[H:58])([H:59])([H:60])[H:61]>>[H][C:3]([C:2]([C:1]([H:59])([H:60])[H:61])([H:57])[H:58])([C:4]([H:54])([H:55])[H:56])[O:5][H]",
        "[C:1]([C@:2]([N:3]([C:4](=[O:5])[O:6][C:7]([C:8]1:[C:9]([H:27]):[C:10]([H:28]):[C:11]([H:29]):[C:12]([H:30]):[C:13]:1[H:31])([H:32])[H:33])[H:34])([C:14](=[O:15])[O:16][C:17]1:[C:18]([H:35]):[C:19]([H:36]):[C:20]([N+:21](=[O:22])[O-:23]):[C:24]([H:37]):[C:25]:1[H:38])[H:39])([H:40])([H:41])[H:42]>>[H]O[C:14]([C@:2]([C:1]([H:40])([H:41])[H:42])([N:3]([C:4](=[O:5])[O:6][C:7]([C:8]1:[C:9]([H:27]):[C:10]([H:28]):[C:11]([H:29]):[C:12]([H:30]):[C:13]:1[H:31])([H:32])[H:33])[H:34])[H:39])=[O:15].[H][O:16][C:17]1:[C:18]([H:35]):[C:19]([H:36]):[C:20]([N+:21](=[O:22])[O-:23]):[C:24]([H:37]):[C:25]:1[H:38]",
        # With At
        "[#6:1](-[#6@:2](-[#8:3]-[At])(-[#6:4](=[#8:5])-[#8:6]-[At:17])-[At])(-[At:19])(-[At:20])-[At:21]>>[#6:1](-[#6:2](=[#8:3])-[#6:4](=[#8:5])-[#8:6]-[At:17])(-[At:19])(-[At:20])-[At:21]",
        "[#8](-[#6](-[#6@@]1(-[At])-[#8]-[#6@](-[#8:7]-[#6@:8]2(-[At:28])-[#6@:9](-[#6:10](-[#8:11]-[At:29])(-[At:30])-[At:31])(-[At:32])-[#8:12]-[#6@:13](-[#8:14]-[At:33])(-[#6:15](-[#8:16]-[At:34])(-[At:35])-[At:36])-[#6@:17]-2(-[#8:18]-[At:37])-[At:38])(-[At])-[#6@@](-[#8]-[At])(-[At])-[#6@](-[#8]-[At])(-[At])-[#6@]-1(-[#8]-[At])-[At])(-[At])-[At])-[At]>>[#8:7](-[#6@:8]1(-[At:28])-[#6@:9](-[#6:10](-[#8:11]-[At:29])(-[At:30])-[At:31])(-[At:32])-[#8:12]-[#6@:13](-[#8:14]-[At:33])(-[#6:15](-[#8:16]-[At:34])(-[At:35])-[At:36])-[#6@:17]-1(-[#8:18]-[At:37])-[At:38])-[At]",
        "CN>>C=O"
    ]

    for rxn in reactions_with_cofactors:
        assert contains_cofactor(rxn) is True, f"Should contain cofactor: {rxn}"

    for rxn in reactions_without_cofactors:
        assert contains_cofactor(rxn) is False, f"Should NOT contain cofactor: {rxn}"
