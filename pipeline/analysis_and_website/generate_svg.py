import os
import csv
from indigo import Indigo, inchi
from indigo.renderer import IndigoRenderer
import time
from pipeline.config import load_paths

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
        pass

def generate_all_svgs(evodex_type, metadata, images_dir):
    if evodex_type in ['F', 'M']:
        return

    full_csv_path = metadata['filename']
    smirks_column = 'smirks'
    try:
        with open(full_csv_path, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                smirks = row[smirks_column]
                evodex_id = row['id']
                svg_filename = f"{evodex_id}-{evodex_type}.svg"
                generate_svg(smirks, svg_filename, images_dir)
    except FileNotFoundError:
        pass

def generate_svgs_for_data_preparation(data_paths, images_dir):
    for key, csv_path in data_paths.items():
        try:
            with open(csv_path, mode='r') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    smirks = row['mapped'] if key == 'raw' else row['smirks']
                    evodex_id = row['rxn_idx'] if key == 'raw' else row['id']
                    svg_filename = f"{evodex_id}-{key}.svg"
                    generate_svg(smirks, svg_filename, images_dir)
        except FileNotFoundError:
            pass

if __name__ == "__main__":
    paths = load_paths('pipeline/config/paths.yaml')
    data_paths = {
        'raw': paths['raw_data'],
        'filtered': paths['filtered_data'],
        'astatine': paths['astatine_data']
    }
    ro_metadata = {
        'R': {'filename': paths['evodex_r']},
        'P': {'filename': paths['evodex_p']},
        'E': {'filename': paths['evodex_e']},
        'N': {'filename': paths['evodex_n']},
        'C': {'filename': paths['evodex_c']},
        'Em': {'filename': paths['evodex_em']},
        'Nm': {'filename': paths['evodex_nm']},
        'Cm': {'filename': paths['evodex_cm']}
    }
    
    generate_svgs_for_data_preparation(data_paths, paths['images_dir'])
    for evodex_type, metadata in ro_metadata.items():
        generate_all_svgs(evodex_type, metadata, paths['images_dir'])
