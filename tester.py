# can be deleted, playground for debugging

from evodex.operators import extract_operator

raw_smiles = '[CH3:1][O:2][C:3](=[O:4])[C@H:5]([CH2:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1)[NH:13][C:14](=[O:15])[C@H:16]([CH2:17][C:18](=[O:19])[OH:20])[NH:21][C:22](=[O:23])[O:24][CH2:25][c:26]1[cH:27][cH:28][cH:29][cH:30][cH:31]1.[OH2:32]>>[C:14](=[O:15])([C@H:16]([CH2:17][C:18](=[O:19])[OH:20])[NH:21][C:22](=[O:23])[O:24][CH2:25][c:26]1[cH:27][cH:28][cH:29][cH:30][cH:31]1)[OH:32].[CH3:1][O:2][C:3](=[O:4])[C@H:5]([CH2:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1)[NH2:13]'

# convert to astatine

# Extract ERO operator (EVODEX-E)
ero_smirks = extract_operator(
    evodex_r_smirks,
    include_stereochemistry=True,
    include_sigma=True,
    include_pi=True,
    include_unmapped_hydrogens=True,
    include_unmapped_heavy_atoms=True,
    include_static_hydrogens=True
)

print("EVODEX-E (ERO) SMIRKS:")
print(ero_smirks)