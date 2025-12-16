# can be deleted, playground for debugging

from evodex.evaluation import operator_matches_reaction

# # Ketone reduction
# rxn = "CC(=O)CCCC>>CC(O)CCCC"
# ro = "[O:45]=[C:46]([C:47])[C:49]>>[H][O:45][C@@:46]([H])([C:47])[C:49]"
# result = operator_matches_reaction(ro, rxn)
# print(result)

# # O-methylation
# rxn = "[H][O:38][C:37]1:[C:35]([O:36][H:65]):[C:34]([H:64]):[C:33]([H:63]):[C:32]([C:31](=[C:30]([C:29]([O:28][H:72])([H:70])[H:71])[H:69])[H:68]):[C:39]:1[H:67]>>[H]C([H])([H])[O:38][C:37]1:[C:35]([O:36][H:65]):[C:34]([H:64]):[C:33]([H:63]):[C:32]([C:31](=[C:30]([C:29]([O:28][H:72])([H:70])[H:71])[H:69])[H:68]):[C:39]:1[H:67]"
# ro = "[H][O:12][C:11]>>[H]C([H])([H])[O:12][C:11]"
# result = operator_matches_reaction(ro, rxn)
# print(result)

ro = "[H][#8:12][#6:11]>>[H]C([H])([H])[#8:12][#6:11]"

print("O-methylation of glycerol, terminal OH")
rxn = "OCC(O)CO>>OCC(O)COC"
result = operator_matches_reaction(ro, rxn)
print(result)

print("O-methylation of glycerol, internal OH")
rxn = "OCC(O)CO>>OCC(OC)CO"
result = operator_matches_reaction(ro, rxn)
print(result)

print("O-methylation of cyclohexanol")
rxn = "C1CCCCC1O>>C1CCCCC1OC"
result = operator_matches_reaction(ro, rxn)
print(result)

print("O-methylation of phenol")
rxn = "C1=CC=CC=C1O>>C1=CC=CC=C1OC"
result = operator_matches_reaction(ro, rxn)
print(result)