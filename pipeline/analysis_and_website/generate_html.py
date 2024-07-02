import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader

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

    # Create index page content
    index_content = []

    # Process each reaction based on EVODEX-R
    evodex_r_df = pd.read_csv(os.path.join(data_dir, 'EVODEX-R_full_reactions.csv'))
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

        context = {
            'evodex_id': evodex_id,
            'reaction_details': reaction_details,
            'smirks': smirks  # Add the original smirks to the context
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
    generate_html_pages(paths, paths['data_dir'], paths['images_dir'], paths['pages_dir'], ['P', 'F', 'M', 'E', 'Em', 'C', 'Cm', 'N', 'Nm'])
