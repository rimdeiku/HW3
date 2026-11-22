# RNR-seq analizė su Snakemake

Homework Part III. Žmogaus RNR-seq duomenų analizė naudojant Snakemake.

## Duomenys

Duomenų rinkinys: GSE52778 (žmogaus kvėpavimo takų raumenų ląstelės, apdorotos deksametazonu).

- SRR1039508 (GSM1275862): untreated
- SRR1039520 (GSM1275874): untreated
- SRR1039521 (GSM1275875): dexamethasone
- SRR1039522 (GSM1275876): albuterol

## Referencija

GENCODE v50, GRCh38 (primary assembly): genomas (`.fa`), anotacija (`.gtf`) ir transkriptomas (`transcripts.fa.gz`).

## Aplankų struktūra

Snakefile            darbo eiga
config/config.yaml   keliai, parametrai, gijos, atmintis
config/samples.csv   mėginių lentelė (sample, condition, fastq_1, fastq_2)
envs/                Conda aplinkos kiekvienai programai
scripts/             duomenų ir referencijos atsisiuntimo skriptai
run_snakemake.sh     paleidimas per SLURM

## Kaip paleisti

1. Aktyvuoti Snakemake aplinką
conda activate snakemake-env

2. Atsisiųsti duomenis ir referenciją
bash scripts/download_reads.sh
bash scripts/download_reference.sh

3. Patikrinti planą (nieko nevykdo)
snakemake -n --use-conda

4. Paleisti per SLURM
sbatch run_snakemake.sh
