from setuptools import setup, find_packages

setup(
    name='evodex',
    version='0.1',
    packages=find_packages(include=['evodex', 'evodex.*']),
    install_requires=[
        'rdkit',  # List your package dependencies here
        'pandas',
        'numpy',
    ],
    include_package_data=True,  # Include additional files specified in MANIFEST.in
    author='J. Christopher Anderson',
    author_email='jcanderson@berkeley.edu',
    description='A project to process enzymatic reactions',
    url='https://github.com/jcaucb/evodex',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
