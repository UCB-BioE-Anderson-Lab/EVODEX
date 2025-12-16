import csv
import sys
import os
import pipeline.web_generation.generate_svg as make_svg
csv.field_size_limit(sys.maxsize)

def count_sources(input_path, output_path):
    entries = []

    with open(input_path, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            source_count = len(row['sources'].split(','))
            entries.append({
                'id': row['id'],
                'source_count': source_count
            })

    # Sort entries by source count descending
    entries.sort(key=lambda x: x['source_count'], reverse=True)

    with open(output_path, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['id', 'source_count'])
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry)

if __name__ == "__main__":
    input_file = 'data/processed/unfiltered_eros.csv'
    output_file = 'data/processed/ero_source_counts.csv'
    count_sources(input_file, output_file)

    # Generate SVGs for entries with only 1 source
    def generate_svgs_for_singletons(input_path, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        with open(input_path, 'r') as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                sources = row['sources'].split(',')
                if len(sources) == 1:
                    smirks = row['smirks']
                    evodex_id = row['id']
                    filename = f"{evodex_id}.svg"
                    make_svg.generate_svg(smirks, filename, output_dir)

    svg_output_dir = 'data/processed/ero_svgs_singletons'
    generate_svgs_for_singletons(input_file, svg_output_dir)

    # Explicitly render a specified SMIRKS string
    # explicit_smirks = "[#6:3]-[#6:4](-[#85:15])(-[#85:16])-[#85].[#16:7]-[#6:8](-[#6:9](-[#6@:10])(-[#85:29])-[#85:30])(-[#85:31])-[#85:32]>>[#6:3]-[#6:4](-[#6:8](-[#16:7])(-[#85:31])-[#85:32])(-[#85:15])-[#85:16].[#6:9](-[#6@:10])(-[#85:29])(-[#85:30])-[#85]"
    # explicit_smirks = "[CH3:1][C:2](=[O:3])[NH:4][C@H:5]1[C@H:6]([O:7][C@H:8]2[C@H:9]([O:10][CH:11]([CH3:12])[C:13](=[O:14])[NH:15][C@@H:16]([CH3:17])[C:18](=[O:19])[NH:20][C@H:21]([CH2:22][CH2:23][C:24](=[O:25])[NH:26][C@@H:27]([CH2:28][CH2:29][CH2:30][C@@H:31]([NH2:32])[C:33](=[O:34])[OH:35])[C:36](=[O:37])[NH:38][C@H:39]([CH3:40])[C:41](=[O:42])[OH:43])[C:44]([NH2:45])=[O:46])[C@@H:47]([NH:48][C:49](=[O:50])[CH2:51][OH:52])[CH:53]([OH:54])[O:55][C@@H:56]2[CH2:57][OH:58])[O:59][C@H:60]([CH2:61][OH:62])[C@@H:63]([OH:64])[C@@H:65]1[OH:66].[CH3:67][S:68][CH2:69][CH2:70][C@@H:71]([NH2:72])[C:73](=[O:74])[OH:75]>>[CH3:1][C:2](=[O:3])[NH:4][C@H:5]1[C@H:6]([O:7][C@H:8]2[C@H:9]([O:10][CH:11]([CH3:12])[C:13](=[O:14])[NH:15][C@@H:16]([CH3:17])[C:18](=[O:19])[NH:20][C@H:21]([CH2:22][CH2:23][C:24](=[O:25])[NH:26][C@@H:27]([CH2:28][CH2:29][CH2:30][C@@H:31]([NH2:32])[C:33](=[O:34])[OH:35])[C:36](=[O:37])[NH:38][C@H:39]([CH2:40][CH2:69][S:68][CH3:67])[C:41](=[O:42])[OH:43])[C:44]([NH2:45])=[O:46])[C@@H:47]([NH:48][C:49](=[O:50])[CH2:51][OH:52])[CH:53]([OH:54])[O:55][C@@H:56]2[CH2:57][OH:58])[O:59][C@H:60]([CH2:61][OH:62])[C@@H:63]([OH:64])[C@@H:65]1[OH:66].[CH3:70][C@@H:71]([NH2:72])[C:73](=[O:74])[OH:75]"
    explicit_smirks = "[OH2:1].[OH:2][C@H:3]1[C@@H:4]([OH:5])[C@@H:6]([OH:7])[O:8][C@@H:9]([O:10][c:11]2[cH:12][cH:13][c:14]3[cH:15][c:16]([Br:17])[cH:18][cH:19][c:20]3[cH:21]2)[C@@H:22]1[OH:23]>>[OH:10][c:11]1[cH:12][cH:13][c:14]2[cH:15][c:16]([Br:17])[cH:18][cH:19][c:20]2[cH:21]1.[OH:1][C@H:6]1[C@H:4]([OH:5])[C@@H:3]([OH:2])[C@@H:22]([OH:23])[C@@H:9]([OH:8])[O:7]1"
    explicit_filename = "glucose_example.svg"
    make_svg.generate_svg(explicit_smirks, explicit_filename, svg_output_dir)