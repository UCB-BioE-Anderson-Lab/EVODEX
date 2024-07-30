
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
The synthesis module provides tools for assigning EVODEX formulas and matching reaction operators. Below is an example usage:

```python
from evodex.synthesis import assign_evodex_F, match_operators

smirks = "CCCO>>CCC=O"
is_valid_formula = assign_evodex_F(smirks)
print(f"{smirks} matches: {is_valid_formula}")

matching_operators = match_operators(smirks, 'C')
print(f"Matching operators for {smirks}: {matching_operators}")
```

For more detailed usage, refer to the [EVODEX Synthesis Demo](https://colab.research.google.com/drive/16liT8RhMCcRzXa_BVdYX7xgbgVAWK4tA).

### Evaluation
The evaluation module provides tools for evaluating reaction operators and synthesis results. Example usage:

```python
from evodex.evaluation import evaluate_reaction_operators, evaluate_synthesis_results, generate_evaluation_report

evaluate_reaction_operators()
evaluate_synthesis_results()
generate_evaluation_report()
```

For more detailed usage, refer to the [EVODEX Evaluation Demo](https://colab.research.google.com/drive/1IvoaXjtnu7ZSvot_1Ovq3g-h5IVCdSn4).

### Mass Spectrometry
The mass spectrometry module provides tools for processing mass spectrometry data. Example usage:

```python
from evodex.mass_spec import process_mass_spec_data

process_mass_spec_data()
```

For more detailed usage, refer to the [EVODEX Mass Spec Demo](https://colab.research.google.com/drive/1CV5HM9lBy-U-J6nLqBlO6Y1WtCFWP8rX).

## Pipeline
The `run_pipeline.py` script is the highest-level runner script for generating EVODEX. This script is not included in the PyPI distribution but can be accessed by cloning the repository. The pipeline processes raw data files and generates the necessary data for EVODEX's functionalities and website.

### Raw Data File
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
    with open("/content/processed_reactions.csv", "wt") as f_out:
        f_out.write(f_in.read())
```

The repository includes a partial version of this file containing only selected reactions. Running the pipeline with this file will reproduce the same results as using the full version.

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
