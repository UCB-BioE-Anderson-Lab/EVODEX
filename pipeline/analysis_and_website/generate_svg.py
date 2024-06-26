import os
import csv
from indigo import Indigo, inchi
from indigo.renderer import IndigoRenderer

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
        indigo.setOption("render-atom-ids-visible", False)  # Do not show atom indices
        indigo.setOption("render-aam-color", "0, 0, 0")  # Set color for atom mapping numbers to black

        reaction.aromatize()
        reaction.correctReactingCenters()

        # Render to file
        renderer.renderToFile(reaction, os.path.join(images_dir, filename))

    except Exception as e:
        print(f"Error rendering reaction image: {e}")

def generate_all_svgs(ro_metadata, data_dir, images_dir):
    for evodex_type, metadata in ro_metadata.items():
        csv_file = metadata['filename']
        smirks_column = 'operator_smirks'
        with open(os.path.join(data_dir, csv_file), mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                smirks = row[smirks_column]
                evodex_id = row['id']
                svg_filename = f"{evodex_id}.svg"
                generate_svg(smirks, svg_filename, images_dir)
