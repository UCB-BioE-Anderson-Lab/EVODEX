# EVODEX Pipeline Documentation

## Overview

This document describes the current data flow and file handling for the EVODEX pipeline, focusing on Phases 3a, 3b, and 3c.  
The key principle is that **all intermediate files live in `data/processed/`** (with clear phase labeling), and **conversion from At → H is performed only once, when needed**.

---

## Phase 3a - EVODEX-E Mining

**Goal:** Extract EVODEX-E operators from EVODEX-P and prune EVODEX-P accordingly.

**Input:**

- `evodex_p_filtered` (At)

**Output:**

- `evodex_e_phase3a_preliminary` (At)
- `evodex_p_phase3a_retained` (At)

**Notes:**

- No At→H conversion happens here.
- Operators with >1 supporting P source are retained.

---

## Phase 3b - EVODEX-E Trimming (Dominance Pruning)

**Goal:** Perform dominance pruning of EVODEX-E, and trim P, F, and R accordingly.

**Input:**

- `evodex_e_phase3a_preliminary` (At)
- `evodex_p_phase3a_retained` (At)
- `evodex_f_filtered` (no SMIRKS, no conversion needed)
- `evodex_r` (At)

**Steps:**

1. Convert required inputs to H:

    - `evodex_e_phase3a_preliminary` → `evodex_e_phase3a_preliminary_H`
    - `evodex_p_phase3a_retained` → `evodex_p_phase3a_retained_H`

    Both written to `data/processed/`.

2. **Copy H versions to `evodex/data/`** to support `match_operators()`:

    - `evodex/data/EVODEX-E_reaction_operators.csv`
    - `evodex/data/EVODEX-P_reaction_operators.csv`

3. Run `match_operators()` to build the dominance table.

4. Perform dominance pruning:  
    - Determine retained operators (no renumbering yet).  
    - Identify retained P sources and prune F and R accordingly.  
    - Write pruned results to `data/processed/`.

**Outputs:** (all still in `data/processed/`, At versions where applicable)

- `evodex_e_phase3b_final` (At)
- `evodex_p_phase3b_final` (At)
- `evodex_f_phase3b_final` (no SMIRKS, no conversion)
- `evodex_r_phase3b_final` (At)

---

## Phase 3c - Publishing to `evodex/data`

**Goal:** Convert final versions to H and publish to `evodex/data/`.

**Input:**

- `evodex_e_phase3b_final` (At)
- `evodex_p_phase3b_final` (At)
- `evodex_f_phase3b_final` (no conversion needed)
- `evodex_r_phase3b_final` (At)

**Steps:**

1. Assign final EVODEX.1 IDs:  
    - `EVODEX-E_reaction_operators.csv`  
    - `EVODEX-P_reaction_operators.csv`  
    - `EVODEX-R_full_reactions.csv`  
    (sources updated to use assigned EVODEX-P and EVODEX-R IDs)

2. Copy `evodex_f_phase3b_final` as:  
    - `EVODEX-F_unique_formulas.csv`  
    (sources updated to use assigned EVODEX-P IDs)

**Result:** Final published versions in `evodex/data/`, used for website and evaluation.

---

## Summary Principles

✅ **All heavy computation happens on At versions.**  
✅ **Conversion to H happens once — controlled and staged.**  
✅ **evodex/data/ is only updated during controlled phases (when needed).**  
✅ **Processed files serve as canonical, reproducible record of all intermediate states.**  
✅ **No smiles-based operations are repeated — only ID mapping and trimming.**
✅ **Final EVODEX.1 ID assignment happens only in Phase 3c.**

---

## TODO

- [ ] Refactor Phase 3b to fully follow this plan.
- [ ] Implement Phase 3c for controlled publishing.
- [ ] Update `paths.yaml` to clarify phase-based file naming.

---
