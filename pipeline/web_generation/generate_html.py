import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from ast import literal_eval

def generate_html_page(template, filename, context, output_dir):
    html_content = template.render(context)
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, filename), 'w') as file:
        file.write(html_content)

def _parse_sources(sources):
    """
    Parse the sources field from the CSV file.

    Parameters:
    sources (str): The sources field as a string.

    Returns:
    list: A list of source strings.
    """
    sources = "" + str(sources)  # Turn the int to a string, or no effect if already a string
    sources = sources.replace('"', '')  # Remove all double quotes
    return sources.split(',')  # Split by commas

def generate_evodex_r_pages(env, evodex_r_df, source_df, pages_dir):
    template = env.get_template('evodex_r_template.html')
    
    for index, row in evodex_r_df.iterrows():
        evodex_id = row['id']
        smirks = row['smirks']
        
        # Parse the sources column using the _parse_sources function
        sources = _parse_sources(row['sources'])

        # Create a copy of the relevant DataFrame slice
        sources_data = source_df[source_df['rxn_idx'].astype(str).isin(sources)].copy()
        sources_data['protein_refs'] = sources_data['protein_refs'].apply(
            lambda refs: ', '.join(literal_eval(refs)) if isinstance(refs, str) and refs.startswith('[') else refs
        )
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

def generate_evodex_p_pages(env, evodex_p_df, evodex_r_df, evodex_f_df, evodex_m_df, evodex_ro_dfs, ro_metadata, pages_dir):
    template = env.get_template('evodex_p_template.html')
    r_to_smirks = {row['id']: row['smirks'] for index, row in evodex_r_df.iterrows()}

    p_to_f_map = {}
    p_to_m_map = {}
    p_to_ro_map = {}

    for df, value_key, id_key in [
        (evodex_f_df, 'formula', 'evodexF'),
        (evodex_m_df, 'mass', 'evodexM')
    ]:
        for index, row in df.iterrows():
            sources = row['sources'].replace('"', '').split(',')
            value = literal_eval(row[value_key]) if value_key == 'formula' else row[value_key]
            if value_key == 'formula' and isinstance(value, dict):
                value = list(value.items())
            for p_id in sources:
                if value_key == 'formula':
                    if p_id not in p_to_f_map:
                        p_to_f_map[p_id] = {'value': value, 'id': row['id']}
                else:
                    if p_id not in p_to_m_map:
                        p_to_m_map[p_id] = {'value': value, 'id': row['id']}

    # Populate RO map
    for evodex_type, df in evodex_ro_dfs.items():
        for index, row in df.iterrows():
            sources = row['sources'].replace('"', '').split(',')
            value = row['smirks']
            for p_id in sources:
                if p_id not in p_to_ro_map:
                    p_to_ro_map[p_id] = {}
                p_to_ro_map[p_id][evodex_type] = {'value': value, 'id': row['id']}

    for index, row in evodex_p_df.iterrows():
        evodex_id = row['id']
        smirks = row.get('smirks', '')
        sources = row['sources'].replace('"', '').split(',')

        f_data = p_to_f_map.get(evodex_id, {'id': '', 'value': []})
        m_data = p_to_m_map.get(evodex_id, {'id': '', 'value': ''})
        ro_data = p_to_ro_map.get(evodex_id, {})

        full_reactions = {r_id: r_to_smirks.get(r_id, '') for r_id in sources}

        context = {
            'evodex_id': evodex_id,
            'smirks': smirks,
            'f_data': f_data,
            'm_data': m_data,
            'ro_data': ro_data,
            'sources': sources,
            'r_to_smirks': full_reactions,
            'ro_metadata': ro_metadata
        }

        generate_html_page(template, f"{evodex_id}.html", context, pages_dir)

def generate_evodex_f_pages(env, evodex_f_df, evodex_p_df, pages_dir):
    template = env.get_template('evodex_f_template.html')
    for index, row in evodex_f_df.iterrows():
        evodex_id = row['id']
        formula = literal_eval(row['formula'])
        sources = row['sources'].replace('"', '').split(',')

        formula_list = list(formula.items()) if isinstance(formula, dict) else []

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

def generate_evodex_ro_pages(env, evodex_ro_dfs, evodex_p_df, ro_metadata, pages_dir):
    template = env.get_template('evodex_ro_template.html')
    for evodex_type, df in evodex_ro_dfs.items():
        for index, row in df.iterrows():
            evodex_id = row['id']
            smirks = row['smirks']
            sources = row['sources'].replace('"', '').split(',')

            partial_reactions = {}
            for p_id in sources:
                if not evodex_p_df.loc[evodex_p_df['id'] == p_id].empty:
                    partial_reaction = evodex_p_df.loc[evodex_p_df['id'] == p_id, 'smirks'].values[0]
                    partial_reactions[p_id] = partial_reaction

            context = {
                'evodex_id': evodex_id,
                'operator_smirks': smirks,
                'sources': sources,
                'partial_reactions': partial_reactions,
                'ro_metadata': ro_metadata,
                'evodex_type': evodex_type
            }

            generate_html_page(template, f"{evodex_id}.html", context, pages_dir)

def generate_type_index_pages(env, evodex_type_dfs, pages_dir):
    template = env.get_template('type_index_template.html')
    for evodex_type, df in evodex_type_dfs.items():
        evodex_ids = df['id'].tolist()

        context = {
            'type_name': f"EVODEX-{evodex_type}",
            'type_subtitle': f"EVODEX-{evodex_type} Index",
            'evodex_ids': evodex_ids
        }

        generate_html_page(template, f"EVODEX-{evodex_type}_index.html", context, pages_dir)

def generate_main_index_page(env, evodex_types, root_dir, ro_metadata):
    template = env.get_template('main_index_template.html')
    context = {
        'evodex_types': evodex_types,
        'ro_metadata': ro_metadata,
        'google_colab_notebooks': [
            {'title': 'EVODEX Synthesis Demo', 'url': 'https://colab.research.google.com/drive/16liT8RhMCcRzXa_BVdYX7xgbgVAWK4tA'},
            {'title': 'EVODEX Evaluation Demo', 'url': 'https://colab.research.google.com/drive/1IvoaXjtnu7ZSvot_1Ovq3g-h5IVCdSn4'},
            {'title': 'EVODEX Mass Spec Demo', 'url': 'https://colab.research.google.com/drive/1CV5HM9lBy-U-J6nLqBlO6Y1WtCFWP8rX'}
        ],
        'additional_links': [
            {'title': 'Correlating EC Number with EVODEX', 'url': 'website/pages/correlating_ec_number.html'},
            {'title': 'Synthesis Subset of EVODEX-E', 'url': 'website/pages/synthesis_subset.html'},
            {'title': 'Mass Spectrometry Subset of EVODEX-M', 'url': 'website/pages/mass_spec_subset.html'}
        ]
    }
    generate_html_page(template, 'index.html', context, root_dir)

def generate_synthesis_subset_page(env, evodex_e_synthesis_df, pages_dir):
    template = env.get_template('synthesis_subset_template.html')
    rows = evodex_e_synthesis_df.to_dict(orient='records')
    for row in rows:
        row['sources'] = [source.strip() for source in row['sources'].split(',')]

    context = {
        'data': rows
    }

    generate_html_page(template, 'synthesis_subset.html', context, pages_dir)

def generate_mass_spec_subset_page(env, evodex_m_subset_df, pages_dir):
    template = env.get_template('mass_spec_subset_template.html')
    rows = evodex_m_subset_df.to_dict(orient='records')
    for row in rows:
        row['sources'] = [source.strip() for source in row['sources'].split(',')]

    context = {
        'data': rows
    }

    generate_html_page(template, 'mass_spec_subset.html', context, pages_dir)

def generate_html_pages(paths, data_dir, images_dir, pages_dir, evodex_types):
    env = Environment(loader=FileSystemLoader(paths['template_dir']))

    evodex_r_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_r'])))
    evodex_p_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_p'])))
    evodex_f_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_f'])))
    evodex_m_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_m'])))
    evodex_e_synthesis_df = pd.read_csv(os.path.join(data_dir, os.path.basename(paths['evodex_e_synthesis'])))
    evodex_m_subset_df = pd.read_csv(os.path.join(data_dir, 'EVODEX-M_mass_spec_subset.csv'))

    evodex_ro_dfs = {}
    evodex_type_dfs = {}

    for evodex_type in evodex_types:
        if evodex_type not in ['R', 'P', 'F', 'M']:
            evodex_ro_dfs[evodex_type] = pd.read_csv(os.path.join(data_dir, os.path.basename(paths[f'evodex_{evodex_type.lower()}'])))
        evodex_type_dfs[evodex_type] = pd.read_csv(os.path.join(data_dir, os.path.basename(paths[f'evodex_{evodex_type.lower()}'])))

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

    source_df = pd.read_csv(paths['raw_data'])

    generate_evodex_r_pages(env, evodex_r_df, source_df, pages_dir)
    generate_evodex_p_pages(env, evodex_p_df, evodex_r_df, evodex_f_df, evodex_m_df, evodex_ro_dfs, ro_metadata, pages_dir)
    generate_evodex_f_pages(env, evodex_f_df, evodex_p_df, pages_dir)
    generate_evodex_m_pages(env, evodex_m_df, evodex_p_df, pages_dir)
    generate_evodex_ro_pages(env, evodex_ro_dfs, evodex_p_df, ro_metadata, pages_dir)
    generate_type_index_pages(env, evodex_type_dfs, pages_dir)
    generate_main_index_page(env, evodex_types, os.path.abspath(os.path.join(pages_dir, os.pardir, os.pardir)), ro_metadata)
    generate_synthesis_subset_page(env, evodex_e_synthesis_df, pages_dir)
    generate_mass_spec_subset_page(env, evodex_m_subset_df, pages_dir)  # Add this line

if __name__ == "__main__":
    paths = {
        'template_dir': 'pipeline/web_generation/templates',
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
        'evodex_cm': 'website/data/EVODEX-Cm_reaction_operators.csv',
        'evodex_e_synthesis': 'website/data/EVODEX-E_synthesis_subset.csv',
        'evodex_m_mass_spec_subset': 'website/data/EVODEX-M_mass_spec_subset.csv',
        
    }
    data_dir = 'website/data'
    images_dir = 'website/images'
    pages_dir = 'website/pages'
    evodex_types = ['R', 'P', 'E', 'N', 'C', 'Em', 'Nm', 'Cm', 'F', 'M']

    generate_html_pages(paths, data_dir, images_dir, pages_dir, evodex_types)
