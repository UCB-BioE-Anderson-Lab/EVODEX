
# EVODEX

EVODEX is a Python package that provides tools for the curation, validation, and data-driven prediction of enzymatic reactions. This package includes core algorithms for synthesis, evaluation, and mass spectrometry analysis, and can be installed via PyPI. Additionally, users can clone the repository to access the full pipeline for data mining and website generation.

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
    - [Synthesis](#synthesis)
    - [Evaluation](#evaluation)
    - [Mass Spectrometry](#mass-spectrometry)
4. [Pipeline](#pipeline)
5. [Data](#data)
6. [Developer Notes](#developer-notes)
7. [Citing EVODEX](#citing-evodex)
8. [License](#license)

## Overview
EVODEX provides two primary modes of use:
1. **PyPI Distribution:** Install EVODEX via pip to access the core algorithms for enzymatic reaction synthesis, evaluation, and mass spectrometry analysis.
2. **Cloning the Repository:** Clone the EVODEX repository to access the complete pipeline for data mining and website generation.

## Installation
To install EVODEX via PyPI, use the following command:
```bash
pip install evodex
```

## Usage

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

## Pipeline
The `run_pipeline.py` script is the highest-level runner script for generating EVODEX. This script is not included in the PyPI distribution but can be accessed by cloning the repository. The pipeline processes raw data files and generates the necessary data for EVODEX's functionalities and website. It includes additional dependencies not required of the PyPI distribution.  See requirements.txt for details. 

### Raw Data File
EVODEX.0, the current distribution of operators, was built from the dataset derived from BRENDA discussed in:
"EnzymeMap: Curation, validation and data-driven prediction of enzymatic reactions" by E. Heid, D. Probst, W. H. Green and G. K. H. Madsen.

To use the full version of the raw data file, download and decompress it as follows:

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
    with open("/content/processed_reactions.csv", "wt") as f_out):
        f_out.write(f_in.read())
```

The repository includes a partial version of this file containing only selected reactions. Running the pipeline with this file will reproduce the same results as using the full version.

## GitHub Pages Website
A partial set of the most common reaction operators (about 100 of each type) generated by the software is available on the GitHub Pages website. This partial dataset provides a preview of the most common operators. The full dataset can be rendered by running `run_pipeline.py`. Access the website here: [EVODEX GitHub Pages](https://ucb-bioe-anderson-lab.github.io/EVODEX/).


## Data
The pipeline initially writes preliminary mined data files and errors to `EVODEX/data`. These files are then reprocessed and saved in the `EVODEX/evodex/data` and `EVODEX/website/data` folders. The files in `evodex/data` are included in the PyPI distribution, while the files in `website/data` are used to generate the EVODEX website.

## Developer Notes
To run some of the main methods in this project in Visual Studio Code, you may need to enable developer mode from the command line. This typically involves setting environment variables or modifying your VS Code configuration to include the necessary paths and dependencies.

## Citing EVODEX
If you use EVODEX, please cite our publication:
"Extraction of Enzymatic Partial Reaction Operators for Biochemical Analysis and Synthesis" by <insert all authors>, and J. Christopher Anderson.

## License
EVODEX is released under the MIT License.

---

### Notebooks
The following Jupyter notebooks demonstrate the usage of the PyPI distribution for the three primary use cases:
- [EVODEX Synthesis Demo](https://colab.research.google.com/drive/16liT8RhMCcRzXa_BVdYX7xgbgVAWK4tA)
- [EVODEX Evaluation Demo](https://colab.research.google.com/drive/1IvoaXjtnu7ZSvot_1Ovq3g-h5IVCdSn4)
- [EVODEX Mass Spec Demo](https://colab.research.google.com/drive/1CV5HM9lBy-U-J6nLqBlO6Y1WtCFWP8rX)

These notebooks are also available in the `notebooks` section of the repository.
