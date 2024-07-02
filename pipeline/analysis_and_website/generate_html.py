import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from pipeline.config import load_paths

def generate_html_pages(paths, data_dir, images_dir, pages_dir, evodex_types):
    os.makedirs(pages_dir, exist_ok=True)
    env = Environment(loader=FileSystemLoader(paths['template_dir']))
    template = env.get_template('reaction_template.html')
    index_template = env.get_template('index_template.html')

    # Load all EVODEX dataframes
    evodex_dataframes = {}
    for evodex_type in evodex_types:
        file_path = os.path.join(data_dir, f"EVODEX-{evodex_type}_reaction_operators.csv")
        if os.path.exists(file_path):
            evodex_dataframes[evodex_type] = pd.read_csv(file_path)

    # Load additional dataframes
    raw_reactions_df = pd.read_csv(paths['raw_data'])
    filtered_reactions_df = pd.read_csv(paths['filtered_data'])
    astatine_reactions_df = pd.read_csv(paths['astatine_data'])
    evodex_r_df = pd.read_csv(paths['evodex_r'])
    evodex_p_df = pd.read_csv(paths['evodex_p'])
    evodex_f_df = pd.read_csv(paths['evodex_f'])
    evodex_m_df = pd.read_csv(paths['evodex_m'])

    # Create index page content
    index_content = []

    # Process each reaction based on EVODEX-R
    for index, row in evodex_r_df.iterrows():
        evodex_id = str(int(row['id']))  # Ensure ID is in string format
        smirks = row['smirks']

        reaction_details = {}
        for evodex_type in evodex_types:
            df = evodex_dataframes.get(evodex_type)
            if df is not None:
                reaction_data = df[df['id'] == int(evodex_id)]  # Ensure matching IDs as integers
                if not reaction_data.empty:
                    reaction_details[evodex_type] = {
                        'smirks': reaction_data.iloc[0]['smirks'],
                        'svg_filename': f"{evodex_id}-{evodex_type}.svg"  # Correct path
                    }

        # Load additional reaction data
        raw_reaction = raw_reactions_df[raw_reactions_df['rxn_idx'] == int(evodex_id)]
        filtered_reaction = filtered_reactions_df[filtered_reactions_df['id'] == int(evodex_id)]
        astatine_reaction = astatine_reactions_df[astatine_reactions_df['id'] == int(evodex_id)]
        formula_reaction = evodex_f_df[evodex_f_df['id'] == int(evodex_id)]
        mass_reaction = evodex_m_df[evodex_m_df['id'] == int(evodex_id)]

        context = {
            'evodex_id': evodex_id,
            'raw_reaction': raw_reaction.iloc[0]['mapped'] if not raw_reaction.empty else 'N/A',
            'raw_svg': f"{evodex_id}-raw.svg",
            'filtered_reaction': filtered_reaction.iloc[0]['smirks'] if not filtered_reaction.empty else 'N/A',
            'filtered_svg': f"{evodex_id}-filtered.svg",
            'astatine_reaction': astatine_reaction.iloc[0]['smirks'] if not astatine_reaction.empty else 'N/A',
            'astatine_svg': f"{evodex_id}-astatine.svg",
            'reaction_details': reaction_details,
            'smirks': smirks,  # Add the original smirks to the context
            'formula': formula_reaction.iloc[0]['smirks'] if not formula_reaction.empty else 'N/A',
            'mass': mass_reaction.iloc[0]['sources'] if not mass_reaction.empty else 'N/A'
        }

        html_content = template.render(context)
        page_filename = f"{evodex_id}.html"
        with open(os.path.join(pages_dir, page_filename), 'w') as file:
            file.write(html_content)

        index_content.append((evodex_id, page_filename))

    # Generate the index page
    index_context = {
        'reactions': index_content
    }
    with open(os.path.join(pages_dir, 'index.html'), 'w') as file:
        file.write(index_template.render(index_context))

if __name__ == "__main__":
    paths = load_paths('pipeline/config/paths.yaml')
    generate_html_pages(paths, paths['data_dir'], paths['images_dir'], paths['pages_dir'], ['P', 'E', 'N', 'C', 'Em', 'Nm', 'Cm'])

