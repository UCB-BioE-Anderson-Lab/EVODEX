

# Reactivity Series Dataset and Pipeline

This directory contains curated quantitative organic reactivity data drawn from primary literature and textbook references. Each dataset is defined in `tables_manifest.yaml` and includes SMIRKS-style reaction templates, rate or equilibrium data, and annotations describing mechanistic context.  The pipeline interprets these tables and generates reaction operators and images for the published analysis.

---

## 1. Overview of the Data

Each entry in `tables_manifest.yaml` specifies:
- **table_id** — short unique slug for the dataset
- **reaction_template** — minimal operator to predict product and compute atom maps
- **reactivity_type** and **units** — quantitative property (rate, pKa, etc.)
- **context** — mechanistic rationale or interpretation
- **source_ref** — reference key linking to the citation list below

| Mechanistic Group | Example Tables |
|--------------------|----------------|
| **SN2** | `steric_sn2`, `resonance_stabilized_sn2`, `alpha_substituent_effects_sn2` |
| **SNAr** | `2chloro5nitropyrimidine_snar` |
| **Electrophilic addition** | `alkene_hydration`, `styrene_hydration`, `alkene_bromination` |
| **Nucleophilic acyl substitution** | `pnpa_amine_amide_formation_a/b`, `trifluoroacetanilide_methanolysis` |
| **Reductions** | `ketone_reduction` |
| **Radical reactions** | `peroxide_radical_abstraction` |
| **Electrophilic aromatic substitution** | `aryl_bromination` |
| **Pericyclic reactions** | `cyclopentadiene_diels_alder` |
| **Acid-base** | `alcohol_pka` |

---

## 2. Processing Workflow

### 2.1 Data Projection
`run_reactivity_series.py`

Reads the YAML tables and generates substrate–product pairs by applying each `reaction_template` to all listed substrates.

**Output:** intermediate JSON with `substrate_smiles`, `product_smiles`, and associated metadata.

---

### 2.2 Deterministic Atom Mapping
`twinmap.py`

Transfers atom-mapping numbers from the operator to the substrate and product:
- Operator atoms define the initial mapping.
- Matched atoms are numbered sequentially.
- Remaining identical atoms (outside the reaction center) receive maps by synchronized BFS propagation.
- Only heavy atoms are assigned (H omitted).

**Output:** mapped substrate/product SMILES and JSON with mapping details.

---

### 2.3 Hydrogen Mapping
`hydrogen_mapping.py`

Converts mapped pairs into EVODEX-P format involving:
- Hydrogen restoration
- Conversion of Hydrogen to Astatine
- Mapping of Astatines
- Conversion back to Hydrogen form

**Output:** Fully-mapped SMIRKS data for operator extraction.

### 2.4 Operator Extraction
`operator_extraction.py`

Calculates the reaction operators for each reaction for:
- EVODEX-A, B, C, D, E
- RdChiral

**Output:** Operator SMIRKS in JSON

### 2.5 Graphics Generation
`graphics_builder.py`

For each table, generates an SVG for inclusion in the analysis writeup with:
- Rendering of the original template reaction operator
- Rendering of each substrate
- Rendering of the substrate portion of each derived operator alongside the substrate
- Inclusion of rate value

**Output:** SVG files

---

## 3. Proposed File and Script Structure

```
tbd
```

---

## 4. References

| Slug | Full Reference |
|------|----------------|
| **Ballinger_1960** | Ballinger, P.; Long, F. A. *J. Am. Chem. Soc.* **1960**, 82 (4), 795–798. DOI 10.1021/ja01489a008. |
| **Jencks_1968** | Jencks, W. P.; Gilchrist, M. *J. Am. Chem. Soc.* **1968**, 90 (10), 2622–2637. DOI 10.1021/ja01012a030. |
| **Campodonico_2020** | Campodónico, P. R.; Olivares, B.; Tapia, R. A. *Front. Chem.* **2020**, 8, 583. DOI 10.3389/fchem.2020.00583. |
| **Howard_Ingold_2011** | Howard, J. A.; Ingold, K. U. *Can. J. Chem.* **1967**, 45 (8), 793–802 [re-cited as Howard & Ingold 2011]. DOI 10.1139/v67-132. |
| **carey_sundberg_volA_3ed** | Carey, F. A.; Sundberg, R. J. *Advanced Organic Chemistry*, Vol. A: *Structure and Mechanisms*, 3rd ed.; Plenum Press: New York, 1990. |
| **carey_sundberg_volB_3ed** | Carey, F. A.; Sundberg, R. J. *Advanced Organic Chemistry*, Vol. B: *Reactions and Synthesis*, 3rd ed.; Plenum Press: New York, 1990. |
| **march_4ed** | March, J. *Advanced Organic Chemistry: Reactions, Mechanisms, and Structure*, 4th ed.; John Wiley & Sons: New York, 1992. |