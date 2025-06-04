# can be deleted, playground for debugging

from evodex.operators import extract_operator

# Example EVODEX-P SMILES
# evodex_p_smirks = "[At:14]-[#7:1]=[#6:2](-[#7:3](-[At:15])-[At:16])-[#7:4](-[At:17])-[#6:5](-[At:18])(-[At:19])-[#6:6](-[At:20])(-[At:21])-[#6:7](-[At:22])(-[At:23])-[#6@:8](-[At:24])(-[#7:9](-[At:25])-[At:26])-[#6:10](=[#8:11])-[#8:12]-[At:27]>>[At:14]-[#7:1](-[At])-[#6:2](-[#7:4](-[At:17])-[#6:5](-[At:18])(-[At:19])-[#6:6](-[At:20])(-[At:21])-[#6:7](-[At:22])(-[At:23])-[#6@:8](-[At:24])(-[#7:9](-[At:25])-[At:26])-[#6:10](=[#8:11])-[#8:12]-[At:27])=[#8].[At:15]-[#7:3](-[At:16])-[At]"
evodex_p_smirks = "[At:19]-[#6:1](-[At:20])(-[At:21])-[#6@:2](-[At:22])(-[#8:3]-[At:23])-[#6:4](=[#8:5])-[#6:6]1=[#7:7]-[#6:8]2:[#6:9](:[#7:10](-[At:24]):[#6:11](=[#7:12]-[At:25]):[#7:13](-[At:26]):[#6:14]:2=[#8:15])-[#7:16](-[At:27])-[#6:17]-1(-[At:28])-[At:29]>>[At:19]-[#6:1](-[At:20])(-[At:21])-[#6@:2](-[At:22])(-[#8:3]-[At:23])-[#6:4](=[#8:5])-[#6:6]1=[#7:7]-[#6:8]2:[#6:9](:[#7:10](-[At:24]):[#6:11](=[#8]):[#7:13](-[At:26]):[#6:14]:2=[#8:15])-[#7:16](-[At:27])-[#6:17]-1(-[At:28])-[At:29].[At:25]-[#7:12](-[At])-[At]"

# Extract ERO operator (EVODEX-E)
ero_smirks = extract_operator(
    evodex_p_smirks,
    include_stereochemistry=True,
    include_sigma=True,
    include_pi=True,
    include_unmapped_hydrogens=True,
    include_unmapped_heavy_atoms=True,
    include_static_hydrogens=True
)

print("EVODEX-E (ERO) SMIRKS:")
print(ero_smirks)