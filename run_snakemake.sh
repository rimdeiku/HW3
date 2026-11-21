#!/bin/bash
#SBATCH --job-name=rnaseq_smk
#SBATCH --partition=main
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

# Naudojimas (iš repozitorijos šaknies, su aktyvuota Snakemake aplinka):
#   conda activate <snakemake aplinka>
#   sbatch run_snakemake.sh

set -euo pipefail

CONFIG="${1:-config/config.yaml}"

cd "${SLURM_SUBMIT_DIR}"
command -v snakemake >/dev/null || { echo "snakemake nerastas: aktyvuokite aplinką prieš sbatch"; exit 1; }

mkdir -p logs
snakemake --use-conda \
          --configfile "${CONFIG}" \
          -c "${SLURM_CPUS_PER_TASK}" \
          "${@:2}"
