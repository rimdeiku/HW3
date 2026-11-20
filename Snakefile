import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"])["sample"].tolist()
OUT = config["outdir"]


rule all:
    input:
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


rule fastqc_raw:
    input:
        f"data/raw/{{sample}}_{{read}}.fastq.gz",
    output:
        html=f"{OUT}/fastqc_raw/{{sample}}_{{read}}_fastqc.html",
        zip=f"{OUT}/fastqc_raw/{{sample}}_{{read}}_fastqc.zip",
    threads: config["threads"]["fastqc"]
    conda:
        "envs/fastqc.yaml"
    log:
        f"logs/fastqc_raw/{{sample}}_{{read}}.log",
    shell:
        "fastqc -t {threads} -o {OUT}/fastqc_raw {input} &> {log}"


rule fastp:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz",
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
        f"logs/fastp/{{sample}}.log",
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
        f"logs/fastqc_trimmed/{{sample}}_{{read}}.log",
    shell:
        "fastqc -t {threads} -o {OUT}/fastqc_trimmed {input} &> {log}"
