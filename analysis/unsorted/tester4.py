import os
from indigo import Indigo
from indigo.renderer import IndigoRenderer
import pandas as pd
from evodex.operators import extract_operator
from evodex.astatine import convert_dataframe_smirks_column_at_to_h

# Setup Indigo
indigo = Indigo()
renderer = IndigoRenderer(indigo)

indigo.setOption("render-output-format", "svg")
indigo.setOption("render-implicit-hydrogens-visible", True)  # Show conventional –OH in carboxylate
indigo.setOption("embedding-uniqueness", "none")
indigo.setOption("render-bond-length", 40)
indigo.setOption("render-atom-ids-visible", False)
indigo.setOption("render-aam-color", "0, 0, 0")
indigo.setOption("render-coloring", True)
indigo.setOption("render-stereo-style", "none")
indigo.setOption("render-label-mode", "hetero")  # Only show labels for heteroatoms

def generate_svg(smirks, filename, images_dir):
    try:
        reaction = indigo.loadReactionSmarts(smirks)

        indigo.setOption("render-output-format", "svg")
        indigo.setOption("render-implicit-hydrogens-visible", False)
        indigo.setOption("embedding-uniqueness", "none")
        indigo.setOption("render-bond-length", 40)
        indigo.setOption("render-atom-ids-visible", False)
        indigo.setOption("render-aam-color", "0, 0, 0")

        reaction.aromatize()
        reaction.correctReactingCenters()

        renderer.renderToFile(reaction, os.path.join(images_dir, filename))
    except Exception as e:
        print(f"[!] Could not render SMIRKS: {e}")

# Define Michael addition SMIRKS
rxn_smiles = "[C:11]([At:16])([At:17])=[C:12]([At:18])[C:13](=[O:14])[O:15]([At:19])>>[C]([At])([At])([At])[S][C:11]([At:16])([At:17])[C:12]([At])([At:18])[C:13](=[O:14])[O:15]([At:19])"

# Define settings for C/N/E abstraction levels
configs = {
    "EVODEX-C": dict(
        include_stereochemistry=True,
        include_sigma=False,
        include_pi=False,
        include_unmapped_hydrogens=True,
        include_unmapped_heavy_atoms=True,
        include_static_hydrogens=True
    ),
    "EVODEX-N": dict(
        include_stereochemistry=True,
        include_sigma=True,
        include_pi=False,
        include_unmapped_hydrogens=True,
        include_unmapped_heavy_atoms=True,
        include_static_hydrogens=True
    ),
    "EVODEX-E": dict(
        include_stereochemistry=True,
        include_sigma=True,
        include_pi=True,
        include_unmapped_hydrogens=True,
        include_unmapped_heavy_atoms=True,
        include_static_hydrogens=True
    )
}

# Extract operators
raw_operators = []
for label, settings in configs.items():
    try:
        smirks = extract_operator(rxn_smiles, **settings)
        raw_operators.append((label, smirks))
    except Exception as e:
        print(f"[!] Error extracting {label}: {e}")

# Convert [At] → [H]
df_raw = pd.DataFrame(raw_operators, columns=["Operator", "SMIRKS"])
df_converted, errors = convert_dataframe_smirks_column_at_to_h(df_raw, column_name="SMIRKS")

# Output directory for SVGs
output_dir = "svg_outputs"
os.makedirs(output_dir, exist_ok=True)

# Generate SVG of the original reaction
try:
    original_smirks_with_h = "[C:11]([H:16])([H:17])=[C:12]([H:18])[C:13](=[O:14])[O:15]([H:19])>>[C]([H])([H])([H])[S][C:11]([H:16])([H:17])[C:12]([H])([H:18])[C:13](=[O:14])[O:15]([H:19])"
    original_rxn = indigo.loadReactionSmarts(original_smirks_with_h)
    indigo.setOption("render-output-format", "svg")
    indigo.setOption("render-implicit-hydrogens-visible", False)
    indigo.setOption("embedding-uniqueness", "none")
    indigo.setOption("render-bond-length", 40)
    indigo.setOption("render-atom-ids-visible", False)
    indigo.setOption("render-aam-color", "0, 0, 0")
    original_rxn.aromatize()
    original_rxn.layout()
    original_rxn.correctReactingCenters()
    renderer.renderToFile(original_rxn, os.path.join(output_dir, "Original_Reaction.svg"))
except Exception as e:
    print(f"[!] Could not render original reaction: {e}")

# Generate SVG for each abstraction level
for _, row in df_converted.iterrows():
    label = row["Operator"]
    smirks = row["SMIRKS"]
    filename = f"{label}.svg"
    generate_svg(smirks, filename, output_dir)

# List of test substrates (name, SMILES)
substrates = [
    ("Butene", 'CC=CC'),
    ("Propene", 'C=CC'),
    ("Propionic_Acid", 'CCC(=O)O'),
    ("Propanol", 'CCCO'),
    ("Butenol", 'C=CCO'),
    ("Acrylate", 'C=CC(=O)O'),
    ("Acrylamide", 'C=CC(=O)NC'),
    ("Acrolein", 'C=CC=O'),
    ("Proprionaldehyde", 'CCC=O'),
    ("methyl propyl ether", 'CCCOC'),
    ("Butane", 'CCCC'),
]

# Generate SVGs for each substrate
for name, smiles in substrates:
    try:
        mol = indigo.loadMolecule(smiles)
        mol.aromatize()
        mol.layout()  # Ensure 2D coordinates are generated
        indigo.setOption("render-implicit-hydrogens-visible", True)
        indigo.setOption("render-label-mode", "hetero")
        renderer.renderToFile(mol, os.path.join(output_dir, f"{name}.svg"))
    except Exception as e:
        print(f"[!] Could not render substrate {name}: {e}")