

import os
from indigo import Indigo
from indigo.renderer import IndigoRenderer

# Setup Indigo
indigo = Indigo()
renderer = IndigoRenderer(indigo)

# Directory to save SVGs
images_dir = "structure_svgs"
os.makedirs(images_dir, exist_ok=True)

# Molecules to render (name, SMILES)
molecules = [
    ("propanol", "CCC[OH]"),
    ("propanal", "CCC=[O]"),
    ("propionic_acid", "CCC(=O)[OH]"),
    ("acrylate", "C=CC(=O)[OH]")
]

def generate_molecule_svg(name, smiles):
    try:
        mol = indigo.loadMolecule(smiles)

        # Set rendering options
        indigo.setOption("render-coloring", False)
        indigo.setOption("render-bond-thickness", 2.0)
        indigo.setOption("render-font-size", 10)
        indigo.setOption("render-label-mode", "hetero")
        indigo.setOption("render-bond-length", 50)  # slightly larger for clarity
        indigo.setOption("render-image-size", 300, 200)  # ensure aspect ratio fits horizontal layout

        mol.aromatize()

        filename = os.path.join(images_dir, f"{name}.svg")
        renderer.renderToFile(mol, filename)
        print(f"Rendered: {filename}")
    except Exception as e:
        print(f"Error rendering {name}: {e}")

# Generate SVGs for all molecules
for name, smiles in molecules:
    generate_molecule_svg(name, smiles)