#!/usr/bin/env bash
# Parsisiunčia GENCODE v50 (GRCh38) genomą, anotaciją ir transkriptomą.
# Paleisti iš repozitorijos šaknies:  sbatch scripts/download_reference.sh
# Reikia: wget.

#SBATCH --job-name=hw3_reference
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.log

set -euo pipefail

GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50"

mkdir -p reference

for file in GRCh38.primary_assembly.genome.fa.gz \
            gencode.v50.primary_assembly.annotation.gtf.gz \
            gencode.v50.transcripts.fa.gz; do
    wget --continue --directory-prefix reference "${GENCODE_URL}/${file}"
done
