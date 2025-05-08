import csv
import os
from collections import defaultdict
from evodex.astatine import hydrogen_to_astatine_reaction
from evodex.mapping import map_atoms
from evodex.utils import reaction_hash
from pipeline0.config import load_paths
from pipeline0.version import __version__


def ensure_directories(paths: dict):
    """Ensure that all necessary directories exist."""
    for path in paths.values():
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

def write_row(writer, data):
    """Write a row to the CSV file."""
    writer.writerow(data)

def process_raw_data(input_file, output_file, ec_representation):
    """Process the initial raw data file."""
    total_raw_reactions = 0
    valid_reactions = 0

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'smirks', 'sources', 'error', 'ec_num']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            total_raw_reactions += 1
            ec_num = row.get('ec_num', '')
            if ec_num:
                ec_representation[ec_num] = False
            try:
                # This is a validation step to exclude any malformed data
                reaction_hash(row['mapped'])
                valid_reactions += 1

                # Write the validated data row
                new_row = {'id': row['rxn_idx'], 'smirks': row['mapped'], 'sources': row.get('rxn_idx', ''), 'error': '', 'ec_num': ec_num}
            except Exception as e:
                # Write an erroneous data row
                new_row = {'id': row['rxn_idx'], 'smirks': '', 'sources': '', 'error': str(e), 'ec_num': ec_num}
            write_row(writer, new_row)

    return total_raw_reactions, valid_reactions

def process_data(input_file, output_file, transformation_function):
    """General function to process data with transformation."""
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['id', 'smirks', 'sources', 'error', 'ec_num']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            if row['error']:  # Skip processing if there's an existing error
                write_row(writer, row)
                continue
            try:
                new_data = transformation_function(row)
                write_row(writer, new_data)
            except Exception as e:
                row['error'] = str(e)
                write_row(writer, row)

def consolidate_reactions(input_file, output_file, raw_input_file, evodex_raw_output_file, ec_representation):
    """Consolidate similar reactions based on reaction_hash and retain up to 4 unique instances of each ec_num."""
    ec_num_map = defaultdict(lambda: defaultdict(dict))
    total_evodex_r_reactions = 0
    evodex_r_ids = set()

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            try:
                smirks = row['smirks']
                ec_num = row.get('ec_num', '')
                if smirks and ec_num and len(ec_num_map[ec_num]) < 4:  # Ensure smirks and ec_num are not empty and limit to 4 instances
                    rxn_hash = reaction_hash(smirks)
                    if rxn_hash not in ec_num_map[ec_num]:  # Ensure unique reaction hash
                        ec_num_map[ec_num][rxn_hash] = row
                        evodex_r_ids.add(row['id'])
            except Exception as e:
                print(f"Error processing row: {row} - {e}")

    sorted_ec_nums = sorted(ec_num_map.keys(), key=lambda ec: list(map(int, ec.split('.'))))

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ['id', 'smirks', 'sources', 'ec_num', 'error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        evodex_id_counter = 1
        for ec_num in sorted_ec_nums:
            for rxn_hash, data in ec_num_map[ec_num].items():
                evodex_id = f'EVODEX.{__version__}-R{evodex_id_counter}'
                data['id'] = evodex_id
                sources = data['sources']
                writer.writerow({
                    'id': evodex_id,
                    'smirks': data['smirks'],
                    'sources': sources,
                    'ec_num': data['ec_num'],
                    'error': data['error']
                })
                ec_representation[ec_num] = True
                evodex_id_counter += 1
                total_evodex_r_reactions += 1
    print(f"Total reactions written: {evodex_id_counter - 1}")

    # Write the filtered raw reactions to evodex_raw_output_file
    with open(raw_input_file, 'r') as infile, open(evodex_raw_output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if row['rxn_idx'] in evodex_r_ids:
                writer.writerow(row)

    return total_evodex_r_reactions

def write_ec_representation(ec_representation, output_file):
    """Write the EC number representation dictionary to a CSV file."""
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['ec_num', 'represented'])
        for ec_num, represented in ec_representation.items():
            writer.writerow([ec_num, represented])

def main():
    paths = load_paths('pipeline/config/paths.yaml')
    ensure_directories(paths)

    ec_representation = {}

    # Process initial raw data
    total_raw_reactions, valid_reactions = process_raw_data(paths['raw_data'], paths['filtered_data'], ec_representation)

    # Subsequent processing steps
    process_data(paths['filtered_data'], paths['astatine_data'], lambda row: {
        'id': row['id'],
        'smirks': hydrogen_to_astatine_reaction(row['smirks']),
        'sources': row['sources'],
        'ec_num': row['ec_num'],
        'error': ''
    })

    process_data(paths['astatine_data'], paths['mapped_data'], lambda row: {
        'id': row['id'],
        'smirks': map_atoms(row['smirks']),
        'sources': row['sources'],
        'ec_num': row['ec_num'],
        'error': ''
    })

    # Consolidate reactions to generate EVODEX-R data with EC number constraints
    total_evodex_r_reactions = consolidate_reactions(
        paths['mapped_data'],
        paths['evodex_r'],
        paths['raw_data'],
        paths['selected_data'],  # Use the selected_data path from paths.yaml
        ec_representation
    )

    # Write EC number representation to a CSV file
    write_ec_representation(ec_representation, os.path.join(paths['errors_dir'], 'ec_representation.csv'))

    # Calculate statistics
    percentage_loss = ((total_raw_reactions - valid_reactions) / total_raw_reactions) * 100
    ec_with_reactions = sum(1 for represented in ec_representation.values() if represented)

    # Print final statistics
    print(f"Total raw reactions processed: {total_raw_reactions}")
    print(f"Valid reactions retained after hash sanitization: {valid_reactions}")
    print(f"Percentage data loss due to hash sanitization: {percentage_loss:.2f}%")
    print(f"Total unique EC numbers encountered: {len(ec_representation)}")
    print(f"EC numbers populated with at least one reaction: {ec_with_reactions}")
    print(f"Total EVODEX-R reactions: {total_evodex_r_reactions}")

if __name__ == "__main__":
    main()
