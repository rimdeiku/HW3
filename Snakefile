import pandas as pd

configfile: "config/config.yaml"

samples_df = pd.read_csv(config["samples"]).set_index("sample", drop=False)
samples = samples_df["sample"].tolist()

wildcard_constraints:
    read="[12]",


def get_fastq(wildcards):
    return samples_df.loc[wildcards.sample, f"fastq_{wildcards.read}"]


def get_fastq_pair(wildcards):
    return {
        "r1": samples_df.loc[wildcards.sample, "fastq_1"],
        "r2": samples_df.loc[wildcards.sample, "fastq_2"],
    }

OUT = config["outdir"]
LOG = config["logdir"]

HISAT2_INDEX = config["indexes"]["hisat2_index"]
HISAT2_FILES = multiext(
    HISAT2_INDEX,
    ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
    ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2",
)

rule all:
    input:
	# 2.1 Sekų nuskaitymų pirminis apdorojimas
        expand(
            f"{OUT}/fastqc_raw/{{sample}}_{{read}}_fastqc.html",
            sample=samples,
            read=[1, 2],
        ),
        expand(
            f"{OUT}/trimmed/{{sample}}_{{read}}.trimmed.fastq.gz",
            sample=samples,
            read=[1, 2],
        ),
        expand(
            f"{OUT}/fastqc_trimmed/{{sample}}_{{read}}.trimmed_fastqc.html",
            sample=samples,
            read=[1, 2],
        ),
	# 2.2 Lygiavimu pagrįstas kiekybinis įvertinimas
        expand(f"{OUT}/aligned/{{sample}}.sorted.bam.bai", sample=samples),
        f"{OUT}/featurecounts/counts.txt",


rule fastqc_raw:
    input:
        get_fastq,
    output:
        html=f"{OUT}/fastqc_raw/{{sample}}_{{read}}_fastqc.html",
        zip=f"{OUT}/fastqc_raw/{{sample}}_{{read}}_fastqc.zip",
    threads: config["threads"]["fastqc"]
    conda:
        "envs/fastqc.yaml"
    log:
        f"{LOG}/fastqc_raw/{{sample}}_{{read}}.log",
    shell:
        """
        tmp=$(mktemp -d)
        ln -s "$(realpath {input})" $tmp/{wildcards.sample}_{wildcards.read}.fastq.gz
        fastqc -t {threads} -o {OUT}/fastqc_raw \
               $tmp/{wildcards.sample}_{wildcards.read}.fastq.gz &> {log}
        rm -r $tmp
        """

rule fastp:
    input:
        unpack(get_fastq_pair),
    output:
        r1=f"{OUT}/trimmed/{{sample}}_1.trimmed.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_2.trimmed.fastq.gz",
        html=f"{OUT}/fastp/{{sample}}.html",
        json=f"{OUT}/fastp/{{sample}}.json",
    params:
        q=config["fastp"]["quality"],
        l=config["fastp"]["length_required"],
    threads: config["threads"]["fastp"]
    conda:
        "envs/fastp.yaml"
    log:
        f"{LOG}/fastp/{{sample}}.log",
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --detect_adapter_for_pe \
              -q {params.q} -l {params.l} -w {threads} \
              -h {output.html} -j {output.json} &> {log}
        """


rule fastqc_trimmed:
    input:
        f"{OUT}/trimmed/{{sample}}_{{read}}.trimmed.fastq.gz",
    output:
        html=f"{OUT}/fastqc_trimmed/{{sample}}_{{read}}.trimmed_fastqc.html",
        zip=f"{OUT}/fastqc_trimmed/{{sample}}_{{read}}.trimmed_fastqc.zip",
    threads: config["threads"]["fastqc"]
    conda:
        "envs/fastqc.yaml"
    log:
        f"{LOG}/fastqc_trimmed/{{sample}}_{{read}}.log",
    shell:
        "fastqc -t {threads} -o {OUT}/fastqc_trimmed {input} &> {log}"


# 2.2 Lygiavimu pagrįstas kiekybinis įvertinimas

# Genomo indeksas: Paleidžia tik jei indekso failų nėra
rule hisat2_index:
    input:
        genome=config["reference"]["genome"],
    output:
        HISAT2_FILES,
    params:
        prefix=HISAT2_INDEX,
    threads: config["threads"]["hisat2_build"]
    resources:
        mem_mb=config["memory_mb"]["hisat2_build"],
    conda:
        "envs/hisat2.yaml"
    log:
        f"{LOG}/hisat2_index/hisat2_build.log",
    shell:
        "hisat2-build -p {threads} {input.genome} {params.prefix} &> {log}"


# 2. Suderinimas su genomu (kiekvienam mėginiui), rezultatas iš karto surikiuotas BAM
rule hisat2_align:
    input:
        r1=f"{OUT}/trimmed/{{sample}}_1.trimmed.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_2.trimmed.fastq.gz",
        index=HISAT2_FILES,
    output:
        bam=f"{OUT}/aligned/{{sample}}.sorted.bam",
    params:
        prefix=HISAT2_INDEX,
        sort_threads=config["threads"]["samtools"],
    threads: config["threads"]["hisat2"]
    resources:
        mem_mb=config["memory_mb"]["hisat2"],
    conda:
        "envs/hisat2.yaml"
    log:
        f"{LOG}/hisat2/{{sample}}.log",
    shell:
        """
        hisat2 -p {threads} -x {params.prefix} \
               -1 {input.r1} -2 {input.r2} 2> {log} \
        | samtools sort -@ {params.sort_threads} -o {output.bam} - 2>> {log}
        """

rule bam_index:
    input:
        f"{OUT}/aligned/{{sample}}.sorted.bam",
    output:
        f"{OUT}/aligned/{{sample}}.sorted.bam.bai",
    threads: config["threads"]["samtools"]
    conda:
        "envs/hisat2.yaml"
    log:
        f"{LOG}/samtools_index/{{sample}}.log",
    shell:
        "samtools index -@ {threads} {input} &> {log}"


# 3. Genų kiekybinis įvertinimas
rule featurecounts:
    input:
        bams=expand(f"{OUT}/aligned/{{sample}}.sorted.bam", sample=samples),
        gtf=config["reference"]["annotation"],
    output:
        counts=f"{OUT}/featurecounts/counts.txt",
        summary=f"{OUT}/featurecounts/counts.txt.summary",
    params:
        strand=config["featurecounts"]["strandness"],
        feature=config["featurecounts"]["feature_type"],
        attribute=config["featurecounts"]["attribute"],
    threads: config["threads"]["featurecounts"]
    conda:
        "envs/subread.yaml"
    log:
        f"{LOG}/featurecounts/featurecounts.log",
    shell:
        """
        featureCounts -T {threads} -p --countReadPairs \
                      -s {params.strand} -t {params.feature} -g {params.attribute} \
                      -a {input.gtf} -o {output.counts} {input.bams} &> {log}
        """
