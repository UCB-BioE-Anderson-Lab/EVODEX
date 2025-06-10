# can be deleted, playground for debugging

from evodex.evaluation import operator_matches_reaction

# # Ketone reduction, with hydrogens
# rxn = "CC(=O)CCCC>>CC(O)CCCC"
# ro = "[O:45]=[C:46]([C:47])[C:49]>>[H][O:45][C@@:46]([H])([C:47])[C:49]"
# result = operator_matches_reaction(ro, rxn)
# print(result)

# # Ketone reduction, with hydrogens, with extra water
# rxn = "CC(=O)CCCC>>CC(O)CCCC.O"
# ro = "[O:45]=[C:46]([C:47])[C:49]>>[H][O:45][C@@:46]([H])([C:47])[C:49]"
# result = operator_matches_reaction(rxn, ro)
# print(result)


# Amide formation
rxn = "CC(=O)O.CN>>CNC(C)=O"
ro = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]"
result = operator_matches_reaction(ro, rxn)
print(result)


# Amide formation
rxn = "NC.CC(=O)O>>CNC(C)=O"
ro = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]"
result = operator_matches_reaction(ro, rxn)
print(result)


# Amide formation
rxn = "[H]O[C:12]([C@@:10]([C:2]([C:1]([H:62])([H:63])[H:64])([C:3]([H:52])([H:53])[H:54])[C:4]([O:5][P:6](=[O:7])([O:8][H:55])[O:9][H:56])([H:57])[H:58])([O:11][H:59])[H:61])=[O:13].[H][N:15]([C:16]([C:17]([C:18](=[O:19])[O:20][H:65])([H:66])[H:67])([H:68])[H:69])[H:70]>>[C:1]([C:2]([C:3]([H:52])([H:53])[H:54])([C:4]([O:5][P:6](=[O:7])([O:8][H:55])[O:9][H:56])([H:57])[H:58])[C@@:10]([O:11][H:59])([C:12](=[O:13])[N:15]([C:16]([C:17]([C:18](=[O:19])[O:20][H:65])([H:66])[H:67])([H:68])[H:69])[H:70])[H:61])([H:62])([H:63])[H:64]"
ro = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]"
result = operator_matches_reaction(ro, rxn)
print(result)



# Amide formation
rxn = "[H]O[C:12]([C@@:10]([C:2]([C:1]([H:62])([H:63])[H:64])([C:3]([H:52])([H:53])[H:54])[C:4]([O:5][P:6](=[O:7])([O:8][H:55])[O:9][H:56])([H:57])[H:58])([O:11][H:59])[H:61])=[O:13].[H][N:15]([C:16]([C:17]([C:18](=[O:19])[O:20][H:65])([H:66])[H:67])([H:68])[H:69])[H:70]>>[C:1]([C:2]([C:3]([H:52])([H:53])[H:54])([C:4]([O:5][P:6](=[O:7])([O:8][H:55])[O:9][H:56])([H:57])[H:58])[C@@:10]([O:11][H:59])([C:12](=[O:13])[N:15]([C:16]([C:17]([C:18](=[O:19])[O:20][H:65])([H:66])[H:67])([H:68])[H:69])[H:70])[H:61])([H:62])([H:63])[H:64]"
ro = "[H][N:13]([C:14])[H:27].[H]O[C:10]([C:8])=[O:11]>>[C:8][C:10](=[O:11])[N:13]([C:14])[H:27]"
result = operator_matches_reaction(ro, rxn)
print(result)