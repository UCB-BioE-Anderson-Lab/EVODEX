import csv
from collections import defaultdict, Counter

# Import hydrogen_to_astatine_reaction and map_atoms
from evodex.astatine import hydrogen_to_astatine_reaction
from evodex.mapping import map_atoms
from evodex.operators import extract_operator
from evodex.utils import reaction_hash

def process_top5_enzyme_reactions(enzyme_reactions, output_file):
    """Convert reactions by applying hydrogen-to-astatine and mapping, then write results."""
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['organism', 'protein_refs', 'protein_db', 'ec_num', 'original_mapped', 'converted_smiles']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for enzyme_key, reactions in enzyme_reactions:
            organism, protein_refs, protein_db, ec_num = enzyme_key
            for row in reactions:
                try:
                    orig_smiles = row['mapped']
                    astatined = hydrogen_to_astatine_reaction(orig_smiles)
                    remapped = map_atoms(astatined)
                    writer.writerow({
                        'organism': organism,
                        'protein_refs': protein_refs,
                        'protein_db': protein_db,
                        'ec_num': ec_num,
                        'original_mapped': orig_smiles,
                        'converted_smiles': remapped
                    })
                except Exception as e:
                    print(f"Error processing reaction {row.get('rxn_idx', '')}: {e}")

def read_reaction_data(file_path):
    enzyme_counts = defaultdict(list)

    with open(file_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            enzyme_key = (
                row['organism'].strip(),
                row['protein_refs'].strip(),
                row['protein_db'].strip(),
                row['ec_num'].strip()
            )
            enzyme_counts[enzyme_key].append(row)

    return enzyme_counts

def extract_cm_nm_em_operators(input_file, output_file):
    """Extract Cm, Nm, Em operators and their hashes for the top 5 enzymes."""
    operator_flags = {
        'Cm': dict(include_stereochemistry=False, include_sigma=False, include_pi=False,
                   include_unmapped_hydrogens=False, include_unmapped_heavy_atoms=False, include_static_hydrogens=False),
        'Nm': dict(include_stereochemistry=False, include_sigma=True, include_pi=False,
                   include_unmapped_hydrogens=False, include_unmapped_heavy_atoms=False, include_static_hydrogens=False),
        'Em': dict(include_stereochemistry=False, include_sigma=True, include_pi=True,
                   include_unmapped_hydrogens=False, include_unmapped_heavy_atoms=False, include_static_hydrogens=False),
    }

    with open(input_file, newline='', encoding='utf-8') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['organism', 'protein_refs', 'protein_db', 'ec_num',
                      'original_mapped', 'converted_smiles',
                      'Cm', 'Cm_hash', 'Nm', 'Nm_hash', 'Em', 'Em_hash']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            output_row = {**row}
            for label, flags in operator_flags.items():
                try:
                    op = extract_operator(row['converted_smiles'], **flags)
                    op_hash = reaction_hash(op)
                except Exception:
                    op = ''
                    op_hash = ''
                output_row[label] = op
                output_row[f"{label}_hash"] = op_hash
            writer.writerow(output_row)

def summarize_operator_hashes_by_enzyme(operator_file):
    """Print a summary of unique operator hashes grouped by enzyme."""
    from collections import defaultdict, Counter

    enzyme_bins = defaultdict(lambda: {'Cm': Counter(), 'Nm': Counter(), 'Em': Counter()})
    row_lookup = {}
    with open(operator_file, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            enzyme_key = (row['organism'], row['protein_refs'], row['protein_db'], row['ec_num'])
            for label in ['Cm', 'Nm', 'Em']:
                h = row[f'{label}_hash']
                smirks = row[label]
                if h and smirks:
                    row_lookup[(enzyme_key, h, label)] = smirks
            if row['Cm_hash']:
                enzyme_bins[enzyme_key]['Cm'][row['Cm_hash']] += 1
            if row['Nm_hash']:
                enzyme_bins[enzyme_key]['Nm'][row['Nm_hash']] += 1
            if row['Em_hash']:
                enzyme_bins[enzyme_key]['Em'][row['Em_hash']] += 1

    print("\nOperator Hash Summary by Enzyme\n" + "="*35)
    for idx, (enzyme_key, hashes) in enumerate(enzyme_bins.items(), 1):
        organism, protein_refs, protein_db, ec_num = enzyme_key
        print(f"\n[{idx}] {organism} — EC {ec_num}")
        if protein_refs or protein_db:
            print(f"     Protein: {protein_refs} @ {protein_db}")
        for label in ['Em', 'Cm', 'Nm']:
            print(f"  {label} Operators ({len(hashes[label])} unique):")
            for h, count in hashes[label].most_common():
                smirks = row_lookup.get((enzyme_key, h, label), '')
                print(f"    • {h} (n={count})\n      {smirks}")

def main():
    file_path = 'data/raw/raw_reactions.csv'
    enzyme_counts = read_reaction_data(file_path)

    # Sort enzymes by number of reactions (descending)
    sorted_enzymes = sorted(enzyme_counts.items(), key=lambda item: len(item[1]), reverse=True)

    print("Top 5 most diverse enzymes by unique EC level 3:\n")
    seen_ec_level3 = set()
    diverse_enzymes = []

    # Each enzyme is defined by a unique combination of organism, EC number, and any provided
    # accession identifiers (protein_refs and protein_db). Entries without accession numbers were exlcuded 
    # as they are not necessarily the same protein.
    # This composite key serves as a proxy
    # for a distinct biochemical entity. We identified the 5 enzymes with the largest number
    # of associated reactions, enforcing uniqueness at EC level 3 (e.g., 1.1.1) to ensure broad
    # mechanistic diversity across enzyme classes.
    
    for enzyme_key, reactions in sorted_enzymes:
        ec_full = enzyme_key[3]
        ec_level3 = ".".join(ec_full.split(".")[:3]) if ec_full else ""
        protein_refs = enzyme_key[1]
        if ec_level3 not in seen_ec_level3 and ec_level3 and protein_refs.strip("[]").strip():
            seen_ec_level3.add(ec_level3)
            diverse_enzymes.append((enzyme_key, reactions))
        if len(diverse_enzymes) == 5:
            break

    for i, (enzyme_key, reactions) in enumerate(diverse_enzymes, 1):
        organism, protein_refs, protein_db, ec_num = enzyme_key
        print(f"{i}. Organism: {organism}, Protein Refs: {protein_refs}, DB: {protein_db}, EC: {ec_num} — {len(reactions)} reactions")

    # After printing the top 5, process and write output file
    output_file = 'data/processed/top5_enzyme_mapped_reactions.csv'
    process_top5_enzyme_reactions(diverse_enzymes, output_file)
    print(f"\nConverted reaction SMILES written to {output_file}")

    extract_cm_nm_em_operators('data/processed/top5_enzyme_mapped_reactions.csv',
                               'data/processed/top5_enzyme_operators.csv')
    print("Operator extraction complete. Results written to data/processed/top5_enzyme_operators.csv")

    summarize_operator_hashes_by_enzyme('data/processed/top5_enzyme_operators.csv')

if __name__ == "__main__":
    main()