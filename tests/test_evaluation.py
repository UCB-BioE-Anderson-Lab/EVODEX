import pytest
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
from evodex.evaluation import operator_matches_reaction

if __name__ == "__main__":
    pytest.main()


# Test for operator_matches_reaction: basic case
def test_operator_matches_reaction_basic_case():
    from evodex.evaluation import operator_matches_reaction

    # Ketone reduction
    rxn = "CC(=O)CCCC>>CC(O)CCCC"
    ro = "[O:45]=[C:46]([C:47])[C:49]>>[H][O:45][C@@:46]([H])([C:47])[C:49]"
    assert operator_matches_reaction(ro, rxn) is True

    # Ketone reduction with extra water
    rxn = "CC(=O)CCCC>>CC(O)CCCC.O"
    ro = "[O:45]=[C:46]([C:47])[C:49]>>[H][O:45][C@@:46]([H])([C:47])[C:49]"
    assert operator_matches_reaction(ro, rxn) is False

    # Ketone reduction, with water added to RO
    rxn = "CC(=O)CCCC>>CC(O)CCCC"
    ro = "[O:45]=[C:46]([C:47])[C:49]>>[H][O:45][C@@:46]([H])([C:47])[C:49].[O]"
    assert operator_matches_reaction(ro, rxn) is False

    # Amide formation
    rxn = "[H]O[C:12]([C@@:10]([C:2]([C:1]([H:62])([H:63])[H:64])([C:3]([H:52])([H:53])[H:54])[C:4]([O:5][P:6](=[O:7])([O:8][H:55])[O:9][H:56])([H:57])[H:58])([O:11][H:59])[H:61])=[O:13].[H][N:15]([C:16]([C:17]([C:18](=[O:19])[O:20][H:65])([H:66])[H:67])([H:68])[H:69])[H:70]>>[C:1]([C:2]([C:3]([H:52])([H:53])[H:54])([C:4]([O:5][P:6](=[O:7])([O:8][H:55])[O:9][H:56])([H:57])[H:58])[C@@:10]([O:11][H:59])([C:12](=[O:13])[N:15]([C:16]([C:17]([C:18](=[O:19])[O:20][H:65])([H:66])[H:67])([H:68])[H:69])[H:70])[H:61])([H:62])([H:63])[H:64]"
    ro = "[H][N:13]([C:14])[H:27].[H]O[C:10]([C:8])=[O:11]>>[C:8][C:10](=[O:11])[N:13]([C:14])[H:27]"
    assert operator_matches_reaction(ro, rxn) is True

    # Amide formation with switched order of substrates
    rxn = "[H][N:15]([C:16]([C:17]([C:18](=[O:19])[O:20][H:65])([H:66])[H:67])([H:68])[H:69])[H:70].[H]O[C:12]([C@@:10]([C:2]([C:1]([H:62])([H:63])[H:64])([C:3]([H:52])([H:53])[H:54])[C:4]([O:5][P:6](=[O:7])([O:8][H:55])[O:9][H:56])([H:57])[H:58])([O:11][H:59])[H:61])=[O:13]>>[C:1]([C:2]([C:3]([H:52])([H:53])[H:54])([C:4]([O:5][P:6](=[O:7])([O:8][H:55])[O:9][H:56])([H:57])[H:58])[C@@:10]([O:11][H:59])([C:12](=[O:13])[N:15]([C:16]([C:17]([C:18](=[O:19])[O:20][H:65])([H:66])[H:67])([H:68])[H:69])[H:70])[H:61])([H:62])([H:63])[H:64]"
    ro = "[H][N:13]([C:14])[H:27].[H]O[C:10]([C:8])=[O:11]>>[C:8][C:10](=[O:11])[N:13]([C:14])[H:27]"
    assert operator_matches_reaction(ro, rxn) is True

    # # Amide formation expecting 2 substrates, but only 1 substrate present
    # rxn = "NCC(=O)O>>NCC(=O)NCC(=O)O"
    # ro = "[H][N:13]([C:14])[H:27].[H]O[C:10]([C:8])=[O:11]>>[C:8][C:10](=[O:11])[N:13]([C:14])[H:27]"
    # assert operator_matches_reaction(ro, rxn) is True

    # Dominance pruning tests
    rxn = "[C:2](=[C:3][C:4]=[O:5])[C:13]([C:12])=[O:14]>>[H][C:13]([C:2]=[C:3][C:4]=[O:5])([C:12])[O:14][H]"
    ro = "[C:3][C:4](=[O:5])[C:6]>>[H][C:4]([C:3])([O:5][H])[C:6]"
    assert operator_matches_reaction(ro, rxn) is True

    rxn = "[H][C@@:12]([C:11])([C:13]([H])([C:14]=[C:15])[H:89])[C:23]>>[C:11][C:12](=[C:13]([C:14]=[C:15])[H:89])[C:23]"
    ro = "[H][C:2]([C:1])([C@@:3]([H])([C:4])[C:25])[H:96]>>[C:1][C:2](=[C:3]([C:4])[C:25])[H:96]"
    assert operator_matches_reaction(ro, rxn) is True

    # O-methylation
    rxn = "[H][O:38][C:37]1:[C:35]([O:36][H:65]):[C:34]([H:64]):[C:33]([H:63]):[C:32]([C:31](=[C:30]([C:29]([O:28][H:72])([H:70])[H:71])[H:69])[H:68]):[C:39]:1[H:67]>>[H]C([H])([H])[O:38][C:37]1:[C:35]([O:36][H:65]):[C:34]([H:64]):[C:33]([H:63]):[C:32]([C:31](=[C:30]([C:29]([O:28][H:72])([H:70])[H:71])[H:69])[H:68]):[C:39]:1[H:67]"
    ro = "[#6:11]-[#8:12]-[H]>>[#6:11]-[#8:12]-[#6](-[H])(-[H])-[H]"
    assert operator_matches_reaction(ro, rxn) is True
