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
```mermaid
flowchart TD
    A(Input: Protein FASTA & Ligand) --> B[BLASTP Search];
    B --> C[Extract Top 10 Hits];
    C --> D[Convert Hits to YAML for Boltz-2];
    D --> E[Boltz-2: Protein-Ligand Structure Prediction];
    E --> F[PLIP: Protein-Ligand Interaction Profiler];
    F --> G(Vizualisation on PyMOL)

    graph TD;
    A-->B;
    A-->C;
    B-->D;
    C-->D;
    C-->E;
    E-->F;
```
