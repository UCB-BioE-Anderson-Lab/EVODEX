import os
import csv
from indigo import Indigo, inchi
from indigo.renderer import IndigoRenderer
import logging
import time

# Setup Indigo
indigo = Indigo()
indigo_inchi = inchi.IndigoInchi(indigo)
renderer = IndigoRenderer(indigo)

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(filename='svg_generation.log', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

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
        start_time = time.time()
        renderer.renderToFile(reaction, os.path.join(images_dir, filename))
        logging.info(f"SVG generated for {filename} in {time.time() - start_time} seconds")

    except Exception as e:
        logging.error(f"Error rendering reaction image for {filename}: {e}")

def generate_all_svgs(evodex_type, metadata, images_dir):
    full_csv_path = metadata['filename']
    smirks_column = 'smirks'
    try:
        with open(full_csv_path, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                smirks = row[smirks_column]
                evodex_id = row['id']
                svg_filename = f"{evodex_id}-{evodex_type}.svg"
                logging.info(f"Processing SVG for ID: {evodex_id}, file: {svg_filename}")
                generate_svg(smirks, svg_filename, images_dir)
    except FileNotFoundError:
        logging.error(f"File not found: {full_csv_path}")

if __name__ == "__main__":
    setup_logging()
    paths = {'data_dir': 'path_to_data_dir', 'images_dir': 'path_to_images_dir'}
    ro_metadata = {'E': {'filename': 'path_to_evodex_e.csv'}}
    generate_all_svgs('E', ro_metadata, paths['images_dir'])
