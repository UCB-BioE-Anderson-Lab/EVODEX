from setuptools import setup, find_packages

setup(
    name="evodex",
    version="2.0.0",
    packages=find_packages(include=["evodex", "evodex.*"]),
    install_requires=[
        "rdkit-pypi",
        "pandas",
        "numpy",
    ],
    include_package_data=True,
    package_data={
        "evodex": [
            "data/EVODEX-B.csv",
            "data/EVODEX-Bm.csv",
            "data/EVODEX-C.csv",
            "data/EVODEX-Cm.csv",
            "data/EVODEX-D.csv",
            "data/EVODEX-Dm.csv",
            "data/EVODEX-E.csv",
            "data/EVODEX-Em.csv",
            "data/EVODEX-F.csv",
            "data/EVODEX-M.csv",
            "data/EVODEX-M_mass_spec_subset.csv",
            "data/EVODEX-P.csv",
            "data/EVODEX-R.csv",
            "data/EVODEX-D_synthesis_subset.csv",
        ],
    },
    author="J. Christopher Anderson",
    author_email="jcanderson@berkeley.edu",
    description="A project to process enzymatic reactions",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jcaucb/evodex",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
