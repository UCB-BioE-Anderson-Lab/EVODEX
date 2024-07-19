from setuptools import setup, find_packages

setup(
    name='evodex',
    version='0.1.19',
    packages=find_packages(include=['evodex', 'evodex.*']),
    install_requires=[
        'rdkit-pypi',
        'pandas',
        'numpy',
    ],
    include_package_data=True,
    package_data={
        'evodex': [
            'data/EVODEX-C_reaction_operators.csv',
            'data/EVODEX-Cm_reaction_operators.csv',
            'data/EVODEX-E_reaction_operators.csv',
            'data/EVODEX-E_synthesis_subset.csv',
            'data/EVODEX-Em_reaction_operators.csv',
            'data/EVODEX-F_unique_formulas.csv',
            'data/EVODEX-M_mass_spec_subset.csv',
            'data/EVODEX-M_unique_masses.csv',
            'data/EVODEX-N_reaction_operators.csv',
            'data/EVODEX-Nm_reaction_operators.csv',
        ],
    },
    author='J. Christopher Anderson',
    author_email='jcanderson@berkeley.edu',
    description='A project to process enzymatic reactions',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/jcaucb/evodex',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
