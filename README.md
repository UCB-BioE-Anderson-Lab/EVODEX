# EVODEX

EVODEX is a Python package that provides tools for the prediction of mechanistically plausible reaction products, validation of reactions, and mass spectrometry interpretation. It can be installed via PyPI for immediate use in Python projects. Alternatively, users can clone the repository to run the full mining pipeline and generate a customized dataset or website.

# Current Release
This is the EVODEX.1 collection. All IDs start with 'EVODEX.1' and are immutable, ensuring they can be externally referenced without collisions or missing references. Future distributions will be numbered EVODEX.2, EVODEX.3, etc., and may not have reverse compatibility with previous EVODEX.0 IDs. For example, EVODEX.1-E2 may not represent the same SMIRKS as EVODEX.0-E2.

## Table of Contents
1. [Installation](#installation)
2. [Usage](#usage)
    - [Synthesis](#synthesis)
    - [Evaluation](#evaluation)
    - [Mass Spectrometry](#mass-spectrometry)
3. [Website and Dataset Access](#website-and-dataset-access)
4. [Building and Running the Mining Pipeline](#building-and-running-the-mining-pipeline)
5. [Citing EVODEX](#citing-evodex)
6. [License](#license)

## Running EVODEX via PyPi

## Installation
To use EVODEX in a python project, use the following command:
```bash
pip install evodex
```

## Usage
The following Jupyter notebooks demonstrate the usage of the PyPI distribution for the three primary use cases:
- [EVODEX Synthesis Demo](https://colab.research.google.com/drive/16liT8RhMCcRzXa_BVdYX7xgbgVAWK4tA)
- [EVODEX Evaluation Demo](https://colab.research.google.com/drive/1IvoaXjtnu7ZSvot_1Ovq3g-h5IVCdSn4)
- [EVODEX Mass Spec Demo](https://colab.research.google.com/drive/1CV5HM9lBy-U-J6nLqBlO6Y1WtCFWP8rX)

These notebooks are also available in the `notebooks` section of the repository.

### Synthesis
The synthesis module provides tools for predicting reaction products using reaction operators. Below is an example usage:

```python
from evodex.synthesis import project_reaction_operator, project_evodex_operator, project_synthesis_operators

# Specify propanol as the substrate as SMILES
substrate = "CCCO"

# Representation of alcohol oxidation as SMIRKS:
smirks = "[H][C:8]([C:7])([O:9][H])[H:19]>>[C:7][C:8](=[O:9])[H:19]"

# Project the oxidation operator on propanol:
result = project_reaction_operator(smirks, substrate)
print("Direct projection: ", result)

# Specify the dehydrogenase reaction by its EVODEX ID:
evodex_id = "EVODEX.0-E2"

# Apply the dehydrogenase operator to propanol
result = project_evodex_operator(evodex_id, substrate)
print("Referenced EVODEX projection: ", result)

# Project All Synthesis Subset EVODEX-E operators on propanol
result = project_synthesis_operators(substrate)
print("All Synthesis Subset projection: ", result)
```

For more detailed usage, refer to the [EVODEX Synthesis Demo](https://colab.research.google.com/drive/16liT8RhMCcRzXa_BVdYX7xgbgVAWK4tA).

### Evaluation
The evaluation module provides tools for evaluating reaction operators and synthesis results. Below is an example usage:

```python
from evodex.evaluation import assign_evodex_F, match_operators

# Define reaction as oxidation of propanol
reaction = "CCCO>>CCC=O"

# Assign EVODEX-F IDs
assign_results = assign_evodex_F(reaction)
print(assign_results)

# Match reaction operators of type 'E' (or C or N)
match_results = match_operators(reaction, 'E')
print(match_results)
```

For more detailed usage, refer to the [EVODEX Evaluation Demo](https://colab.research.google.com/drive/1IvoaXjtnu7ZSvot_1Ovq3g-h5IVCdSn4).

### Mass Spectrometry
The mass spectrometry module provides tools for predicting masses and identifying reaction operators. Below is an example usage:

```python
from evodex.mass_spec import calculate_mass, find_evodex_m, get_reaction_operators, predict_products

# Calculate exact mass of the compound cortisol as an [M+H]+ ion
cortisol_M_plus_H = "O=C4\C=C2/[C@]([C@H]1[C@@H](O)C[C@@]3([C@@](O)(C(=O)CO)CC[C@H]3[C@@H]1CC2)C)(C)CC4.[H+]"
mass = calculate_mass(cortisol_M_plus_H)

# Define observed masses
substrate_mass = 363.2166 # The expected mass for cortisol's ion
potential_product_mass = 377.2323 # A mass of unknown identity
mass_diff = potential_product_mass - substrate_mass

# Find matching EVODEX-M entries
matching_evodex_m = find_evodex_m(mass_diff, 0.01)
print(matching_evodex_m)

# Get reaction operators
matching_operators = get_reaction_operators(mass_diff, 0.01)
print(matching_operators)

# Predict product structures
predicted_products = predict_products(cortisol_M_plus_H, mass_diff, 0.01)
print(predicted_products)
```

For more detailed usage, refer to the [EVODEX Mass Spec Demo](https://colab.research.google.com/drive/1CV5HM9lBy-U-J6nLqBlO6Y1WtCFWP8rX).

## Website and Dataset Access

A static website for exploring the EVODEX.1 dataset is available here:

🔗 https://ucb-bioe-anderson-lab.github.io/evodex-1-site/

This site provides a browsable, hyperlinked index of all EVODEX.1 operators and includes links to Colab demos and operator definitions.

If you prefer to work offline or programmatically, you can download the full set of operator tables (in CSV format) directly from this directory:

📂 https://github.com/UCB-BioE-Anderson-Lab/EVODEX/tree/main/evodex/data

## Building and Running the Mining Pipeline
The mining pipeline allows you to reproduce the full EVODEX operator set from raw data. This is the process we use to generate the operators that ship with the PyPI distribution. Running the mining pipeline is only required if you want to modify the data, change the curation process, or experiment with new reactions.

⚠️ This project depends on a narrow window of compatible package versions. It is strongly recommended to use Python 3.11.13 and avoid newer Python versions unless you are familiar with rebuilding the RDKit dependency stack using Conda.

### Build Environment

- This release is tested and works with Python **3.11.13**.
- Using other Python versions (especially ≥ 3.12) may cause installation issues with key dependencies like `rdkit-pypi`, `scipy`, and `pandas`.
- Do not attempt to upgrade individual packages or change the Python version unless you're using Conda, which is not recommended for this build.
- For consistency and stability, we recommend using this exact Python version and installing via `venv`.

```bash
/opt/homebrew/bin/python3.11 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## Running the Pipeline

You can run the entire pipeline automatically by running:

```bash
python run_pipeline.py
```

Alternatively, you can run the pipeline step-by-step. First, download and prepare the raw data file, then run the modules sequentially.

To download and prepare the raw data file:

```python
import requests
import gzip

# Download the file
url = "https://github.com/hesther/enzymemap/blob/main/data/processed_reactions.csv.gz?raw=true"
r = requests.get(url)
with open("/content/processed_reactions.csv.gz", "wb") as f:
    f.write(r.content)

# Decompress the file
with gzip.open("/content/processed_reactions.csv.gz", "rt") as f_in:
    with open("/content/processed_reactions.csv", "wt") as f_out:
        f_out.write(f_in.read())
```

Then run the following modules sequentially:

```bash
python -m pipeline.phase1_data_preparation
python -m pipeline.phase2_formula_pruning
python -m pipeline.phase3_ero_mining
python -m pipeline.phase3a_ero_pruning
python -m pipeline.phase3b_ero_trimming
python -m pipeline.phase3c_ero_publishing
python -m pipeline.phase4_operator_completion
python -m pipeline.phase5_mass_subset
python -m pipeline.phase6_synthesis_subset
python -m pipeline.phase7_website
```

## Citing EVODEX
If you use EVODEX, please cite our preprint:
https://www.biorxiv.org/content/10.1101/2025.06.19.660615v1.full

## License
EVODEX is released under the MIT License.
