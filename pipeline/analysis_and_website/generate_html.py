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

    # Function to get reaction data and handle missing columns
    def get_reaction_data(df, column_name, idx_column_name='id'):
        if idx_column_name not in df.columns:
            return 'N/A'
        reaction_data = df[df[idx_column_name] == int(evodex_id)]
        if not reaction_data.empty:
            error_data = reaction_data.iloc[0].get('error', '')
            if pd.notna(error_data) and error_data != '':
                return error_data
            data_value = reaction_data.iloc[0].get(column_name, 'N/A')
            return data_value if pd.notna(data_value) else 'N/A'
        return 'N/A'

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
                if 'id' in df.columns:
                    reaction_data = df[df['id'] == int(evodex_id)]  # Ensure matching IDs as integers
                    if not reaction_data.empty:
                        smirks_data = reaction_data.iloc[0]['smirks']
                        error_data = reaction_data.iloc[0].get('error', '')
                        reaction_details[evodex_type] = {
                            'smirks': error_data if pd.notna(error_data) and error_data != '' else smirks_data,
                            'svg_filename': f"{evodex_id}-{evodex_type}.svg" if not (pd.notna(error_data) and error_data != '') else None
                        }

        # Load additional reaction data
        raw_reaction = get_reaction_data(raw_reactions_df, 'mapped', 'rxn_idx')
        raw_svg = f"{evodex_id}-raw.svg"

        filtered_reaction = get_reaction_data(filtered_reactions_df, 'smirks')
        filtered_svg = f"{evodex_id}-filtered.svg" if filtered_reaction != 'N/A' else None

        astatine_reaction = get_reaction_data(astatine_reactions_df, 'smirks')
        astatine_svg = f"{evodex_id}-astatine.svg" if astatine_reaction != 'N/A' else None

        formula_reaction = get_reaction_data(evodex_f_df, 'formula')
        mass_reaction = get_reaction_data(evodex_m_df, 'mass')

        context = {
            'evodex_id': evodex_id,
            'raw_reaction': raw_reaction,
            'raw_svg': raw_svg,
            'filtered_reaction': filtered_reaction,
            'filtered_svg': filtered_svg,
            'astatine_reaction': astatine_reaction,
            'astatine_svg': astatine_svg,
            'reaction_details': reaction_details,
            'smirks': smirks,  # Add the original smirks to the context
            'formula': formula_reaction,
            'mass': mass_reaction
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
    with open(os.path.join(pages_dir, 'curation_index.html'), 'w') as file:
        file.write(index_template.render(index_context))

if __name__ == "__main__":
    paths = load_paths('pipeline/config/paths.yaml')
    generate_html_pages(paths, paths['data_dir'], paths['images_dir'], paths['pages_dir'], ['R', 'P', 'E', 'N', 'C', 'Em', 'Nm', 'Cm'])
