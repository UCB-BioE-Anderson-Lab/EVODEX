import os
import math
import pandas as pd
import svgwrite
from xml.etree import ElementTree as ET
from indigo import Indigo
from indigo.renderer import IndigoRenderer
from evodex.synthesis import project_synthesis_operators
from rdkit import Chem

# Initialize Indigo
indigo = Indigo()
renderer = IndigoRenderer(indigo)

# Indigo rendering options for SVG export
indigo.setOption("render-output-format", "svg")
indigo.setOption("render-implicit-hydrogens-visible", True)
indigo.setOption("render-label-mode", "hetero")
indigo.setOption("render-coloring", False)
indigo.setOption("render-bond-length", 40)

# Set output directory
output_dir = "figure4_svg"
os.makedirs(output_dir, exist_ok=True)

# Define substrate
substrate = "CCCO"  # Propanol

# Run synthesis projection
results = project_synthesis_operators(substrate)

# Flatten to a list of (evodex_id, product_smiles)
df_result = pd.DataFrame(list(results.items()), columns=["EVODEX ID", "Predicted Products"])
df_result = df_result.explode("Predicted Products")

# Phase 1: Generate and save SVG images
for idx, row in df_result.iterrows():
    evodex_id = row["EVODEX ID"]
    product = row["Predicted Products"]
    try:
        mol = indigo.loadMolecule(product)
        mol.aromatize()
        mol.foldHydrogens()
        mol.layout()
        renderer.renderToFile(mol, os.path.join(output_dir, f"{evodex_id}.svg"))
    except Exception as e:
        print(f"[!] Could not render {evodex_id}: {e}")