#!/usr/bin/env bash
# Parsisiunčia RNA-seq mėginius iš NCBI SRA ir suspaudžia FASTQ.
# Paleisti iš repozitorijos šaknies:  sbatch scripts/download_reads.sh
# Reikia: sra-tools (fasterq-dump), gzip.

#SBATCH --job-name=hw3_reads
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.log

set -euo pipefail

SAMPLES="SRR1039520 SRR1039508 SRR1039521 SRR1039522"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p data/raw

for sample in ${SAMPLES}; do
    if [[ -s "data/raw/${sample}_2.fastq.gz" ]]; then
        continue
    fi
    fasterq-dump "${sample}" --split-3 --force --threads "${THREADS}" --outdir data/raw
    gzip --force data/raw/"${sample}"_*.fastq
done
