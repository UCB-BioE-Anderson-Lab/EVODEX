from setuptools import setup, find_packages

setup(
    name='evodex',
    version='0.1.5',
    packages=find_packages(include=['evodex', 'evodex.*']),
    install_requires=[
        'rdkit-pypi',
        'pandas',
        'numpy',
    ],
    include_package_data=True,
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
