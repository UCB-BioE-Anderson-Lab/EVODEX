import os
from jinja2 import Template
import pandas as pd

def generate_html_pages(paths, data_dir, images_dir, pages_dir):
    template_path = os.path.join(paths['template_dir'], 'reaction_template.html')
    with open(template_path) as template_file:
        template_content = template_file.read()

    template = Template(template_content)

    evodex_p_df = pd.read_csv(paths['evodex_p'])
    r_to_p_map = {}
    for index, row in evodex_p_df.iterrows():
        for r_id in eval(row.get('sources', '[]')):  # Use .get to avoid KeyError
            if r_id not in r_to_p_map:
                r_to_p_map[r_id] = []
            r_to_p_map[r_id].append({
                'id': row['id'],
                'smirks': row['partial_reaction']
            })

    source_filename = paths['raw_data']
    evodex_r_df = pd.read_csv(paths['evodex_r'])

    # Debugging statements to print column names and sample row
    print(f"Columns in {paths['evodex_r']}: {evodex_r_df.columns.tolist()}")
    print(f"Sample row from {paths['evodex_r']}: {evodex_r_df.iloc[0].to_dict()}")

    source_df = pd.read_csv(source_filename)

    for index, row in evodex_r_df.iterrows():
        try:
            evodex_id = row['id']
            smirks = row['mapped_reaction']
            rxn_ids = eval(row.get('sources', '[]'))  # Use .get to avoid KeyError

            sources_data = source_df[source_df['rxn_idx'].isin(rxn_ids)]
            sources = sources_data[['rxn_idx', 'orig_rxn_text', 'natural', 'organism', 'protein_refs', 'protein_db', 'ec_num']].to_dict(orient='records')
            partial_reactions = r_to_p_map.get(evodex_id, [])

            context = {
                'evodex_id': evodex_id,
                'svg_filename': f"{evodex_id}.svg",
                'smirks': smirks,
                'sources': sources,
                'partial_reactions': partial_reactions
            }

            html_content = template.render(context)
            with open(os.path.join(pages_dir, f"{evodex_id}.html"), 'w') as file:
                file.write(html_content)
        except KeyError as e:
            print(f"KeyError: {e} in row: {row.to_dict()}")
