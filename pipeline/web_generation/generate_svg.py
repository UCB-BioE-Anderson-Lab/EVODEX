import os
import csv
from indigo import Indigo, inchi
from indigo.renderer import IndigoRenderer

# Setup Indigo
indigo = Indigo()
indigo_inchi = inchi.IndigoInchi(indigo)
renderer = IndigoRenderer(indigo)

def generate_svg(smirks, filename, images_dir):
    try:
        # Load reaction SMIRKS
        reaction = indigo.loadReactionSmarts(smirks)

        # Set rendering options
        indigo.setOption("render-output-format", "svg")
        indigo.setOption("render-implicit-hydrogens-visible", False)
        indigo.setOption("embedding-uniqueness", "none")
        indigo.setOption("render-bond-length", 40)
        indigo.setOption("render-atom-ids-visible", False)
        indigo.setOption("render-aam-color", "0, 0, 0")

        reaction.aromatize()
        reaction.correctReactingCenters()

        # Render to file
        renderer.renderToFile(reaction, os.path.join(images_dir, filename))

    except Exception as e:
        print(f"Error rendering reaction image: {e}")

def generate_all_svgs(evodex_type, csv_path, images_dir):
    if evodex_type in ['F', 'M']:
        return

    smirks_column = 'smirks'
    
    try:
        with open(csv_path, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                smirks = row[smirks_column]
                evodex_id = row['id']
                svg_filename = f"{evodex_id}.svg"
                generate_svg(smirks, svg_filename, images_dir)
    except FileNotFoundError:
        print(f"File not found: {csv_path}")
    except KeyError:
        print(f"Column '{smirks_column}' not found in the file {csv_path}")
