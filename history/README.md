# EVODEX History

This folder collects earlier versions of the EVODEX project. They document how the operator extraction workflow evolved and how the same ideas were implemented in different chemoinformatics toolkits.

EVODEX began at UC Berkeley with early experiments on enzyme promiscuity and reaction abstractions. The first full implementation was written at 20n using ChemAxon and a combined BRENDA–MetaCyc–KEGG dataset. When the work returned to Berkeley, the code was cleaned up for open-source release. After the academic ChemAxon license expired, the pipeline was reimplemented using the open-source Indigo toolkit. The current EVODEX system is a complete Python rewrite using RDKit.

These historical versions are kept for reference and for anyone interested in adapting operator extraction to other toolkits. They are not maintained, and some depend on proprietary components.

## What’s Here

- **chemaxon/**  
  Last Java version that used ChemAxon. This is the original end-to-end operator extraction pipeline from the 20n era. It requires commercial libraries and is provided for historical context only.

- **indigo/**  
  Java implementation using the open-source Indigo toolkit. This replaced the ChemAxon dependency and reflects the second major generation of the algorithm.

- **History of the EVODEX Project.docx**  
  A full narrative history of the scientific and software development path that led to EVODEX.