# =================================================================== #
# Snakefile
# Description :
#   Snakemake workflow for protein homologous sequence searching using BLAST
#   and structure prediction using Boltz-2.
#
# Author : Samantha SAMSON
# Date : 2026/01/29
# Usage :
#   ~$ bash bash_worflow.sh
# =================================================================== #

#!/usr/bin/env bash
set -euo pipefail  # Si erreur, sortir immédiatement

FASTA_FILE="./data/P00330.fasta"
LIGAND_FILE="./data/NAD_smiles.txt"
HITS_FASTA=""
MAX_TARGETS=10
BLASTP_DIR="./blastp_results"
HITS_DIR="./hits_yalm"
BOLTZ_DIR="./boltz_results"

# Vérification des fichiers d'entrée
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Erreur: Fichier FASTA introuvable: $FASTA_FILE"
    exit 1
fi

if [[ ! -f "$LIGAND_FILE" ]]; then
    echo "Erreur: Fichier ligand introuvable: $LIGAND_FILE"
    exit 1
fi

# 0. Setup des répertoires de résultats.
mkdir -p "$BLASTP_DIR" "$HITS_DIR" "$BOLTZ_DIR"

echo "=== Étape 1: BLASTP ==="
# Lance BLASTP bases de données distantes de NCBI et extrait les 10 premiers hits.
blastp -query "$FASTA_FILE" \
    -db nr -remote \
    -max_target_seqs "$MAX_TARGETS" \
    -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sseq" \
    -evalue 1e-5 \
    -out "$BLASTP_DIR/blastp_results.tsv"

# Vérification des résultats BLAST
if [[ ! -s "$BLASTP_DIR/blastp_results.tsv" ]]; then
    echo "Erreur: Résultats BLASTP vides"
    exit 1
else
    HITS_FASTA="$BLASTP_DIR/hit_sequences.fasta"
    awk 'BEGIN{FS="\t"; OFS="\n"} NR>0 {gsub(/-/, "", $13); if($13!="") print ">"$2, $13}' "$BLASTP_DIR/blastp_results.tsv" > "$HITS_FASTA"
    echo "Séquences hits extraites dans: $HITS_FASTA"
fi

# Vérification des séquences hits
if [[ ! -s "$HITS_FASTA" ]]; then
    echo "Erreur: Aucune séquence extraite du BLAST"
    exit 1
fi

echo -e "\n=== Étape 2: Conversion en YAML pour Boltz-2 ==="
# Convertir les séquences FASTA et le ligand en format YAML pour Boltz-2.
python ./scripts/fasta_to_boltz_yaml.py "$HITS_FASTA" -lig "$LIGAND_FILE" -o "$HITS_DIR"

echo -e "\n=== Étape 3: Prédiction de structure avec Boltz-2 ==="
# Lancer Boltz-2 pour la prédiction des complexes protéine-ligand.
boltz predict "./$HITS_DIR/" --use_msa_server --use_potentials --output_format pdb --accelerator cpu --out_dir "./$BOLTZ_DIR"
# Meilleures predictions avec les options : --recycling_steps 10 --diffusion_samples 25 --> /!\ plus long à exécuter
# Erreur : ValueError: CCD component ASP not found!
# Solution : retélécharger mols.tar sur https://huggingface.co/boltz-community/boltz-2/blob/main/mols.tar et le décompresser dans le dossier cache .boltz

echo -e "\n=== Pipeline terminé! ==="
echo "Résultats dans: $BOLTZ_DIR/predictions/"
