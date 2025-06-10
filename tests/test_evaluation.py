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

    # Ketone reduction, with hydrogens
    rxn = "CC(=O)CCCC>>CC(O)CCCC"
    ro = "[O:45]=[C:46]([C:47])[C:49]>>[H][O:45][C@@:46]([H])([C:47])[C:49]"
    assert operator_matches_reaction(ro, rxn) is True

    # Ketone reduction, with hydrogens, with extra water
    rxn = "CC(=O)CCCC>>CC(O)CCCC.O"
    ro = "[O:45]=[C:46]([C:47])[C:49]>>[H][O:45][C@@:46]([H])([C:47])[C:49]"
    assert operator_matches_reaction(ro, rxn) is False

    # Ketone reduction, with oxygen added to RO
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

    # Amide formation expecting 2 substrates, but only 1 substrate present
    rxn = "NCC(=O)O>>NCC(=O)NCC(=O)O"
    ro = "[H][N:13]([C:14])[H:27].[H]O[C:10]([C:8])=[O:11]>>[C:8][C:10](=[O:11])[N:13]([C:14])[H:27]"
    assert operator_matches_reaction(ro, rxn) is True