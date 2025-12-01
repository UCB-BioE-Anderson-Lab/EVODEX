# Reactivity Series Analysis Using EVODEX Abstractions

This directory contains the code and data used to compare EVODEX abstractions (A–E) and RdChiral operators across a unified collection of reactivity series. The goal is to evaluate how these abstraction levels capture reactivity patterns across spontaneous reactions representing the major mechanistic classes of chemistry. Each table is adapted from experimentally measured reactivity data reported in advanced organic chemistry textbooks or primary literature. The pipeline projects operator templates onto substrates, generates EVODEX-P representations, computes all operator abstractions, and produces the graphics used in the analysis.

---

## File Structure

```
.
├── readme.md
├── tables_manifest.yaml
├── data/
│   └── (reactivity tables in CSV format)
├── 1_project.py
├── 2_ops.py
├── 3_graphics.py
└── projectionmap.py
```

### `tables_manifest.yaml`
Each entry in this manifest defines a single reactivity dataset. For every table, the manifest records:

- a short identifier (`table_id`)
- a descriptive summary of the reaction series
- the operator template to be projected (`reaction_template`, SMIRKS)
- the type of reactivity measurement (`reactivity_type`) and its units
- mechanistic context explaining the observed trend (`context`)
- literature provenance (`source_ref`, `source_page`, `provenance_type`)

These fields allow the pipeline to both execute the projection workflow and retain scientific metadata about the origin and interpretation of each reactivity table.

### `data/`
Contains the reactivity series, one table per mechanism class, with:
- substrate label  
- substrate SMILES  
- experimental reactivity values  

---

## Components of the Pipeline

### 1. `1_project.py` — Template Projection and EVODEX‑P Generation

This script is the entry point. For each table in the manifest it:

- loads the operator template and substrates  
- calls `projectionmap.py` to apply the operator and perform mapping 
- writes the resulting EVODEX‑P representation to `out/`  

The mapped substrate–product pair is treated as the EVODEX‑P abstraction implied by the combination of template and substrate.

---

### 2. `projectionmap.py` — Projection and Mapping Algorithm

This module performs the operator‑projection and mapping logic used by `1_project.py`. Given an unmapped substrate SMILES and an operator SMIRKS, it:

- converts the substrate to an astatine‑labeled form to track hydrogens  
- assigns sequential isotopes across all atoms  
- applies the operator with RDKit and selects the first product  
- promotes isotope labels to atom map numbers in both substrate and product  
- removes map numbers that appear only on one side  
- converts astatine atoms back to hydrogens  
- returns explicit‑hydrogen mapped substrate and product SMILES  

The returned pair defines the EVODEX‑P abstraction for the projected reaction.

---

### 3. `2_ops.py` — Operator Extraction (EVODEX A–E and RdChiral)

Using the EVODEX‑P outputs, this script computes:

- **EVODEX‑A**  
- **EVODEX‑B**  
- **EVODEX‑C**  
- **EVODEX‑D**  
- **EVODEX‑E**  
- **RdChiral** (radius=1, special_groups=True)

All operator outputs are written to `out/`.

---

### 4. `3_graphics.py` — Reaction Series Graphics

This script produces one SVG per reactivity table. Each graphic displays:

- an example full reaction  
- the original substrate  
- the measured reactivity  
- the substrate‑side of the extracted operators

These graphics are the ones included in the supplemental materials.

---

## Running the Pipeline

Run the scripts in order:

1. `python 1_project.py`  
2. `python 2_ops.py`  
3. `python 3_graphics.py`

All output files are written to the `out/` directory, which is created automatically if absent.

---

## References Cited

- **Ballinger_1960**  
  Ballinger, P.; Long, F. A. *Acid Ionization Constants of Alcohols. II. Acidities of Some Substituted Methanols and Related Compounds.* **J. Am. Chem. Soc.** 1960, 82, 795–798. https://doi.org/10.1021/ja01489a008.

- **Campodonico_2020**  
  Campodónico, P. R.; Olivares, B.; Tapia, R. A. *Experimental Analyses Emphasize the Stability of the Meisenheimer Complex in a SNAr Reaction Toward Trends in Reaction Pathways.* **Front. Chem.** 2020, 8, 583. https://doi.org/10.3389/fchem.2020.00583.

- **deGrip_and_Lugtenburg_2022**  
  de Grip, W. J.; Lugtenburg, J. *Isorhodopsin: An Undervalued Visual Pigment Analog.* **Colorants** 2022, 1 (3), 256–279. https://doi.org/10.3390/colorants1030016.

- **Howard_Ingold_2011**  
  Howard, J. A.; Ingold, K. U. *Absolute Rate Constants for Hydrocarbon Autoxidation. VI. Alkyl Aromatic and Olefinic Hydrocarbons.* **Can. J. Chem.** 1967, 45, 793–802. https://doi.org/10.1139/v67-132.

- **Jencks_1968**  
  Jencks, W. P.; Gilchrist, M. *Nonlinear Structure–Reactivity Correlations. The Reactivity of Nucleophilic Reagents toward Esters.* **J. Am. Chem. Soc.** 1968, 90, 2622–2637. https://doi.org/10.1021/ja01012a030.

- **carey_sundberg_volA_3ed**  
  Carey, F. A.; Sundberg, R. J. *Advanced Organic Chemistry,* 3rd ed.; Plenum Press: New York, 1990.

- **carey_sundberg_volB_3ed**  
  Carey, F. A.; Sundberg, R. J. *Advanced Organic Chemistry,* 3rd ed.; Plenum Press: New York, 1990.

- **march_4ed**  
  March, J. *Advanced Organic Chemistry: Reactions, Mechanisms, and Structure,* 4th ed.; Wiley-Interscience: New York, 1992.