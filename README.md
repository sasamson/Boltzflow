# Boltzflow

## Description
**Pipeline automatisé** pour identifier et prédire la structure de complexes protéine-ligand.

Le pipeline suit ces étapes principales :  
1. BLAST pour identifier les séquences homologues.  
2. Prédiction des structures protéine-ligand avec Boltz-2.  
3. Identification des interactions protéine-ligand avec PLIP.  
4. Possibilités d’amélioration pour un pipeline plus robuste et performant.

---

## Pipeline Workflow

+----------------+ +---------------------+ +--------------------------+ +------------------------+
| Input FASTA | ----> | BLAST Search | ----> | Convert hits to YAML | ----> | Boltz-2 Prediction |
| & Ligand(s) | | Extract Top 10 Hits| | (protein + ligand) | | Protein-Ligand Complex |
+----------------+ +---------------------+ +--------------------------+ +------------------------+
|
v
+---------------------------------------------+
| PLIP: Protein-Ligand Interaction Identifier |
+---------------------------------------------+
|
v
+-------------------------------+
|Results visualization on PyMOL |
+-------------------------------+

