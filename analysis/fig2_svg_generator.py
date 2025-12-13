# fig2_svg_generator.py

import os
from indigo import Indigo
from indigo.renderer import IndigoRenderer

SMIRKS = "[C:1]([H:4])([H:5])([H:6])[C:2]([H:7])([H])[O:3][H]>>[C:1]([H:4])([H:5])([H:6])[C:2]([H:7])=[O:3]"
OUT_DIR = "structure_svgs"
OUT_SVG = os.path.join(OUT_DIR, "fig2_operator.svg")

os.makedirs(OUT_DIR, exist_ok=True)

indigo = Indigo()
renderer = IndigoRenderer(indigo)

# Match the options from your Colab snippet
indigo.setOption("render-output-format", "svg")
indigo.setOption("render-implicit-hydrogens-visible", False)
indigo.setOption("embedding-uniqueness", "none")
indigo.setOption("render-bond-length", 40)
indigo.setOption("render-atom-ids-visible", False)
indigo.setOption("render-aam-color", "0, 0, 0")

# Load + prep + render
reaction = indigo.loadReactionSmarts(SMIRKS)
reaction.aromatize()
reaction.correctReactingCenters()

renderer.renderToFile(reaction, OUT_SVG)
print(f"Rendered: {OUT_SVG}")