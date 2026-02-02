# =================================================================== #
# Snakefile
# Description :
#   Snakemake workflow for protein homologous sequence searching using BLAST
#   and structure prediction using Boltz-2.
#
# Author : Samantha SAMSON
# Date : 2026/01/29
# Usage :
#   snakemake --use-conda --cores 4
#   snakemake --use-conda --cores 4 --conda-frontend mamba
# =================================================================== #

import os
from pathlib import Path

# =================================================================== #
# Configuration
# =================================================================== #
# Fichiers d'entrée
FASTA_FILE = "data/P00330.fasta"
LIGAND_FILE = "data/NAD_smiles.txt"

# Paramètres BLAST
MAX_TARGETS = 10
BLAST_EVALUE = 1e-5

# Répertoires de sortie
BLASTP_DIR = "blastp_results"
HITS_DIR = "hits_yaml"
BOLTZ_DIR = "boltz_results_hits_yaml"

# Paramètres Boltz-2
BOLTZ_PARAMS = {
    "use_msa_server": True,
    "use_potentials": True,
    "output_format": "pdb",
    "accelerator": "cpu",  # Changer en "gpu" si disponible
    # Pour meilleures prédictions (plus lent):
    # "recycling_steps": 10,
    # "diffusion_samples": 25
}


def get_yaml_patterns(wildcards):
    """
    Récupère les patterns YAML : 'protid_ligid'
    """
    ckpt = checkpoints.fasta_to_yaml.get()
    yaml_files = Path(ckpt.output.yaml_dir).glob("*.yaml")
    patterns = [f.stem for f in yaml_files]
    return patterns

# =================================================================== #
# Règle finale
# =================================================================== #
rule all:
    input:
        lambda wildcards: [f"{BOLTZ_DIR}/{p}/.done" for p in get_yaml_patterns(wildcards)]

# =================================================================== #
# Étape 0 : Vérification des fichiers d'entrée
# =================================================================== #
rule check_inputs:
    """
    Vérifie que les fichiers d'entrée existent
    """
    input:
        fasta = FASTA_FILE,
        ligand = LIGAND_FILE
    output:
        flag = temp(".inputs_checked")
    run:
        if not Path(input.fasta).exists():
            raise FileNotFoundError(f"Fichier FASTA introuvable: {input.fasta}")
        if not Path(input.ligand).exists():
            raise FileNotFoundError(f"Fichier ligand introuvable: {input.ligand}")
        
        # Vérifier que le FASTA n'est pas vide
        if Path(input.fasta).stat().st_size == 0:
            raise ValueError(f"Fichier FASTA vide: {input.fasta}")
        
        # Créer le flag de validation
        Path(output.flag).touch()
        print("Fichiers d'entrée validés")

# =================================================================== #
# Étape 1 : BLAST search
# =================================================================== #
rule blastp_search:
    """
    Lance BLASTP sur les bases de données distantes de NCBI
    et extrait les 10 premiers hits
    """
    input:
        fasta = FASTA_FILE,
        flag = ".inputs_checked"
    output:
        tsv = f"{BLASTP_DIR}/blastp_results.tsv"
    params:
        max_targets = MAX_TARGETS,
        evalue = BLAST_EVALUE
    threads: 1  # BLAST remote utilise 1 thread
    conda:
        "envs/boltz.yaml"
    shell:
        """
        echo "Query: {input.fasta}"
        echo "Max targets: {params.max_targets}"
        echo "E-value: {params.evalue}"
        
        mkdir {BLASTP_DIR}        
        blastp -query {input.fasta} \
            -db nr -remote \
            -max_target_seqs {params.max_targets} \
            -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sseq" \
            -evalue {params.evalue} \
            -out {output.tsv} \
            || { echo "BLASTP failed"; exit 1; }

        echo "BLAST terminé avec succès !"
        echo "Nombre de hits: $(wc -l < {output.tsv})"
        """

# =================================================================== #
# Étape 2 : Extraction des séquences hits
# =================================================================== #
rule extract_hit_sequences:
    """
    Extrait les séquences des hits BLAST et crée un fichier FASTA
    """
    input:
        tsv = f"{BLASTP_DIR}/blastp_results.tsv"
    output:
        fasta = f"{BLASTP_DIR}/hit_sequences.fasta"
    conda:
        "envs/boltz.yaml"
    shell:
        """
        echo "Extraction des séquences hits..."
        
        awk 'BEGIN{{FS="\\t"; OFS="\\n"}} NR>0 {{
            gsub(/-/, "", $13);
            if($13!="") print ">"$2, $13
        }}' {input.tsv} > {output.fasta}

        echo "Séquences extraites: {output.fasta}"
        echo "Nombre de séquences: $(grep -c "^>" {output.fasta})"
        """

# =================================================================== #
# Étape 3 : Conversion FASTA en YAML pour Boltz-2
# =================================================================== #
checkpoint fasta_to_yaml:
    """
    Convertit les séquences FASTA et le ligand en format YAML pour Boltz-2
    """
    input:
        fasta = f"{BLASTP_DIR}/hit_sequences.fasta",
        ligand = LIGAND_FILE
    output:
        yaml_dir = directory(HITS_DIR)
    conda:
        "envs/boltz.yaml"  # Utilise Biopython
    shell:
        """
        mkdir {HITS_DIR}
        python ./scripts/fasta_to_boltz_yaml.py \
            {input.fasta} \
            -lig {input.ligand} \
            -o {output.yaml_dir}
        
        echo "Conversion terminée: {output.yaml_dir}"
        """

# =================================================================== #
# Étape 4 : Prédiction de structure avec Boltz-2
# =================================================================== #
rule boltz_predict:
    """
    Lance Boltz-2 pour la prédiction des complexes protéine-ligand
    """
    input:
        yaml = lambda wildcards: f"{HITS_DIR}/{wildcards.pattern}.yaml"
    output:
        flag = f"{BOLTZ_DIR}/{{pattern}}/.done"
    params:
        use_msa = "--use_msa_server" if BOLTZ_PARAMS.get("use_msa_server") else "",
        use_potentials = "--use_potentials" if BOLTZ_PARAMS.get("use_potentials") else "",
        output_format = f"--output_format {BOLTZ_PARAMS.get('output_format', 'pdb')}",
        accelerator = f"--accelerator {BOLTZ_PARAMS.get('accelerator', 'cpu')}",
        recycling = f"--recycling_steps {BOLTZ_PARAMS['recycling_steps']}" if 'recycling_steps' in BOLTZ_PARAMS else "",
        diffusion = f"--diffusion_samples {BOLTZ_PARAMS['diffusion_samples']}" if 'diffusion_samples' in BOLTZ_PARAMS else ""
    threads: 8
    resources:
        mem_mb = 16000,
        gpu = 1 if BOLTZ_PARAMS.get("accelerator") == "gpu" else 0
    conda:
        "envs/boltz.yaml"
    shell:
        """
        mkdir -p {BOLTZ_DIR}/{wildcards.pattern}

        # Afficher les paramètres
        echo "Input: {input.yaml}"
        echo "Output: {BOLTZ_DIR}"
        echo "Accelerator: {params.accelerator}"
                
        # Lancer Boltz-2
        boltz predict {input.yaml} \
            {params.use_msa} \
            {params.use_potentials} \
            {params.output_format} \
            {params.accelerator} \
            {params.recycling} \
            {params.diffusion} \
            --out_dir {BOLTZ_DIR}/{wildcards.pattern}

        touch {output.flag}    
        echo "Prédictions terminées"
        echo "Résultats dans: {BOLTZ_DIR}"
        """
