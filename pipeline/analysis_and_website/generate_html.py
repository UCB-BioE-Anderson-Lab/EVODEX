import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from ast import literal_eval

def generate_html_page(template, filename, context, pages_dir):
    html_content = template.render(context)
    os.makedirs(pages_dir, exist_ok=True)
    with open(os.path.join(pages_dir, filename), 'w') as file:
        file.write(html_content)

def generate_evodex_r_pages(env, evodex_r_df, source_df, pages_dir):
    template = env.get_template('evodex_r_template.html')
    for index, row in evodex_r_df.iterrows():
        evodex_id = row['id']
        smirks = row['smirks']
        sources = row['sources'].replace('"', '').split(',')
        sources = [source.strip() for source in sources]

        sources_data = source_df[source_df['rxn_idx'].astype(str).isin(sources)]
        sources_dict = sources_data[['rxn_idx', 'orig_rxn_text', 'natural', 'organism', 'protein_refs', 'protein_db', 'ec_num']].to_dict(orient='records')
        partial_reactions = r_to_p_map.get(evodex_id, [])

        context = {
            'evodex_id': evodex_id,
            'svg_filename': f"{evodex_id}.svg",
            'smirks': smirks,
            'sources': sources_dict,
            'partial_reactions': partial_reactions
        }

        generate_html_page(template, f"{evodex_id}.html", context, pages_dir)

def generate_evodex_p_pages(env, evodex_p_df, evodex_r_df, evodex_f_df, evodex_m_df, additional_evodex_dfs, pages_dir, ro_metadata):
    template = env.get_template('evodex_p_template.html')
    r_to_smirks = {row['id']: row['smirks'] for index, row in evodex_r_df.iterrows()}
    for index, row in evodex_p_df.iterrows():
        evodex_id = row['id']
        smirks = row['smirks']
        sources = row['sources'].replace('"', '').split(',')
        sources = [source.strip() for source in sources]

        f_data = next((item for item in evodex_f_df.itertuples() if evodex_id in item.sources.replace('"', '').split(',')), None)
        m_data = next((item for item in evodex_m_df.itertuples() if evodex_id in item.sources.replace('"', '').split(',')), None)

        ro_data = {}
        for evodex_type in ['E', 'Em', 'N', 'Nm', 'C', 'Cm']:
            df = additional_evodex_dfs[evodex_type]
            ro_data[evodex_type] = [
                {
                    'id': item.id,
                    'smirks': item.smirks
                }
                for item in df.itertuples() if evodex_id in item.sources.replace('"', '').split(',')
            ]

        full_reactions = {r_id: r_to_smirks.get(r_id, '') for r_id in sources}

        context = {
            'evodex_id': evodex_id,
            'svg_filename': f"{evodex_id}.svg",
            'smirks': smirks,
            'sources': sources,
            'f_data': f_data,
            'm_data': m_data,
            'ro_data': ro_data,
            'r_to_smirks': full_reactions,
            'ro_metadata': ro_metadata  # Ensure ro_metadata is passed to the context
        }

        generate_html_page(template, f"{evodex_id}.html", context, pages_dir)

def generate_evodex_f_pages(env, evodex_f_df, evodex_p_df, pages_dir):
    template = env.get_template('evodex_f_template.html')
    for index, row in evodex_f_df.iterrows():
        evodex_id = row['id']
        formula = literal_eval(row['formula'])
        sources = row['sources'].replace('"', '').split(',')
        sources = [source.strip() for source in sources]

        # Ensure formula is a list of tuples for Jinja2 to iterate
        formula_list = list(formula.items()) if isinstance(formula, dict) else []

        # Create a mapping from source ID to details from evodex_p_df
        source_details = []
        for p_id in sources:
            if not evodex_p_df.loc[evodex_p_df['id'] == p_id].empty:
                partial_reaction = evodex_p_df.loc[evodex_p_df['id'] == p_id, 'smirks'].values[0]
                source_details.append({'id': p_id, 'smirks': partial_reaction})

        context = {
            'evodex_id': evodex_id,
            'formula': formula_list,
            'sources': source_details
        }

        generate_html_page(template, f"{evodex_id}.html", context, pages_dir)

def generate_evodex_m_pages(env, evodex_m_df, evodex_p_df, pages_dir):
    template = env.get_template('evodex_m_template.html')
    for index, row in evodex_m_df.iterrows():
        evodex_id = row['id']
        mass = row['mass']
        sources = row['sources'].replace('"', '').split(',')
        sources = [source.strip() for source in sources]

        # Create a mapping from source ID to details from evodex_p_df
        source_details = []
        for p_id in sources:
            if not evodex_p_df.loc[evodex_p_df['id'] == p_id].empty:
                partial_reaction = evodex_p_df.loc[evodex_p_df['id'] == p_id, 'smirks'].values[0]
                source_details.append({'id': p_id, 'smirks': partial_reaction})

        context = {
            'evodex_id': evodex_id,
            'mass': mass,
            'sources': source_details
        }

        generate_html_page(template, f"{evodex_id}.html", context, pages_dir)

def generate_evodex_ro_pages(env, additional_evodex_dfs, evodex_p_df, pages_dir):
    template = env.get_template('evodex_ro_template.html')
    for evodex_type, df in additional_evodex_dfs.items():
        if evodex_type in ['F', 'M']:
            continue
        for index, row in df.iterrows():
            evodex_id = row['id']
            smirks = row['smirks']
            sources = row['sources'].replace('"', '').split(',')
            sources = [source.strip() for source in sources]

            # Create a mapping from source ID to details from evodex_p_df
            partial_reactions = {}
            for p_id in sources:
                if not evodex_p_df.loc[evodex_p_df['id'] == p_id].empty:
                    partial_reaction = evodex_p_df.loc[evodex_p_df['id'] == p_id, 'smirks'].values[0]
                    partial_reactions[p_id] = partial_reaction

            context = {
                'evodex_id': evodex_id,
                'smirks': smirks,
                'sources': sources,
                'partial_reactions': partial_reactions
            }

            generate_html_page(template, f"{evodex_id}.html", context, pages_dir)

def generate_type_index_pages(env, additional_evodex_dfs, pages_dir):
    template = env.get_template('type_index_template.html')
    for evodex_type, df in additional_evodex_dfs.items():
        evodex_ids = df['id'].tolist()

        context = {
            'type_name': f"EVODEX-{evodex_type}",
            'type_subtitle': f"EVODEX-{evodex_type} Index",
            'evodex_ids': evodex_ids
        }

        generate_html_page(template, f"EVODEX-{evodex_type}_index.html", context, pages_dir)

def generate_main_index_page(env, evodex_types, pages_dir, ro_metadata):
    template = env.get_template('main_index_template.html')
    context = {
        'evodex_types': evodex_types,
        'ro_metadata': ro_metadata
    }
    generate_html_page(template, 'index.html', context, pages_dir)

def generate_html_pages(paths, data_dir, images_dir, pages_dir, evodex_types):
    env = Environment(loader=FileSystemLoader(paths['template_dir']))

    # Load DataFrames and ensure rxn_idx is a string
    evodex_r_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_r'])))
    evodex_p_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_p'])))
    evodex_f_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_f'])))
    evodex_m_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_m'])))

    additional_evodex_dfs = {}
    for evodex_type in evodex_types:
        additional_evodex_dfs[evodex_type] = pd.read_csv(os.path.join(data_dir, os.path.basename(paths[f'evodex_{evodex_type.lower()}'])))

    ro_metadata = {
        'E': {'title': 'Reaction Operator E'},
        'Em': {'title': 'Reaction Operator Em'},
        'N': {'title': 'Reaction Operator N'},
        'Nm': {'title': 'Reaction Operator Nm'},
        'C': {'title': 'Reaction Operator C'},
        'Cm': {'title': 'Reaction Operator Cm'}
    }

    global r_to_p_map
    r_to_p_map = {}
    for index, row in evodex_p_df.iterrows():
        sources = row['sources'].replace('"', '').split(',')
        for r_id in sources:
            r_id = r_id.strip()
            if r_id not in r_to_p_map:
                r_to_p_map[r_id] = []
            r_to_p_map[r_id].append({
                'id': row['id'],
                'smirks': row['smirks']
            })

    # Load raw data for EVODEX-R pages
    source_df = pd.read_csv(paths['raw_data'])

    # Generate all types of pages
    generate_evodex_r_pages(env, evodex_r_df, source_df, pages_dir)
    generate_evodex_p_pages(env, evodex_p_df, evodex_r_df, evodex_f_df, evodex_m_df, additional_evodex_dfs, pages_dir, ro_metadata)
    generate_evodex_f_pages(env, evodex_f_df, evodex_p_df, pages_dir)
    generate_evodex_m_pages(env, evodex_m_df, evodex_p_df, pages_dir)
    generate_evodex_ro_pages(env, additional_evodex_dfs, evodex_p_df, pages_dir)
    generate_type_index_pages(env, additional_evodex_dfs, pages_dir)
    generate_main_index_page(env, evodex_types, pages_dir, ro_metadata)

if __name__ == "__main__":
    paths = {
        'template_dir': 'pipeline/analysis_and_website/templates',
        'raw_data': 'website/data/raw_reactions.csv',
        'evodex_r': 'website/data/EVODEX-R_full_reactions.csv',
        'evodex_p': 'website/data/EVODEX-P_partial_reactions.csv',
        'evodex_f': 'website/data/EVODEX-F_unique_formulas.csv',
        'evodex_m': 'website/data/EVODEX-M_unique_masses.csv',
        'evodex_e': 'website/data/EVODEX-E_reaction_operators.csv',
        'evodex_n': 'website/data/EVODEX-N_reaction_operators.csv',
        'evodex_c': 'website/data/EVODEX-C_reaction_operators.csv',
        'evodex_em': 'website/data/EVODEX-Em_reaction_operators.csv',
        'evodex_nm': 'website/data/EVODEX-Nm_reaction_operators.csv',
        'evodex_cm': 'website/data/EVODEX-Cm_reaction_operators.csv'
    }
    data_dir = 'website/data'
    images_dir = 'website/images'
    pages_dir = 'website/pages'
    evodex_types = ['R', 'P', 'E', 'N', 'C', 'Em', 'Nm', 'Cm', 'F', 'M']

    generate_html_pages(paths, data_dir, images_dir, pages_dir, evodex_types)
