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

GENE_BED = config["reference"]["gene_bed"]
GENE_BODY_BED = config["reference"]["gene_body_bed"]

BAMS = expand(f"{OUT}/aligned/{{sample}}.sorted.bam", sample=samples)
BAIS = expand(f"{OUT}/aligned/{{sample}}.sorted.bam.bai", sample=samples)

KALLISTO_INDEX = config["indexes"]["kallisto_index"]

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
        # 2.3 Lygiavimo nereikalaujantis kiekybinis įvertinimas
        expand(f"{OUT}/kallisto/{{sample}}/abundance.tsv", sample=samples),
        # 2.4 Lygiavimo kokybės vertinimas
        expand(f"{OUT}/rseqc/{{sample}}.bam_stat.txt", sample=samples),
        expand(f"{OUT}/rseqc/{{sample}}.infer_experiment.txt", sample=samples),
        expand(f"{OUT}/rseqc/{{sample}}.read_distribution.txt", sample=samples),
        expand(f"{OUT}/rseqc/{{sample}}.geneBodyCoverage.txt", sample=samples),
        expand(f"{OUT}/rseqc/{{sample}}.inner_distance_freq.txt", sample=samples),
        f"{OUT}/deeptools/correlation_heatmap.png",
        f"{OUT}/deeptools/pca.png",


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


# 2.3 Lygiavimo nereikalaujantis kiekybinis įvertinimas

# Kallisto indeksas: paleidžiamas tik jei indekso failo nėra
rule kallisto_index:
    input:
        transcriptome=config["reference"]["transcriptome"],
    output:
        KALLISTO_INDEX,
    threads: config["threads"]["kallisto_index"]
    resources:
        mem_mb=config["memory_mb"]["kallisto_index"],
    conda:
        "envs/kallisto.yaml"
    log:
        f"{LOG}/kallisto_index/kallisto_index.log",
    shell:
        "kallisto index -i {output} {input.transcriptome} &> {log}"


# Transkriptų gausos įvertinimas kiekvienam mėginiui
rule kallisto_quant:
    input:
        r1=f"{OUT}/trimmed/{{sample}}_1.trimmed.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_2.trimmed.fastq.gz",
        index=KALLISTO_INDEX,
    output:
        tsv=f"{OUT}/kallisto/{{sample}}/abundance.tsv",
        info=f"{OUT}/kallisto/{{sample}}/run_info.json",
    params:
        outdir=f"{OUT}/kallisto/{{sample}}",
    threads: config["threads"]["kallisto"]
    resources:
        mem_mb=config["memory_mb"]["kallisto"],
    conda:
        "envs/kallisto.yaml"
    log:
        f"{LOG}/kallisto/{{sample}}.log",
    shell:
        """
        kallisto quant -i {input.index} -o {params.outdir} \
                       -t {threads} {input.r1} {input.r2} &> {log}
        """


# 2.4 Lygiavimo kokybės vertinimas

# BED12 anotacija RSeQC įrankiams (iš GENCODE GTF)
rule gtf_to_bed:
    input:
        gtf=config["reference"]["annotation"],
    output:
        bed=GENE_BED,
    params:
        genepred=GENE_BED + ".genepred",
    conda:
        "envs/ucsc.yaml"
    log:
        f"{LOG}/reference/gtf_to_bed.log",
    shell:
        """
        gtfToGenePred -genePredExt -ignoreGroupsWithoutExons \
                      {input.gtf} {params.genepred} 2> {log}
        genePredToBed {params.genepred} {output.bed} 2>> {log}
        rm -f {params.genepred}
        """

# a. Suderinimo statistika
rule rseqc_bam_stat:
    input:
        bam=f"{OUT}/aligned/{{sample}}.sorted.bam",
    output:
        f"{OUT}/rseqc/{{sample}}.bam_stat.txt",
    resources:
        mem_mb=config["memory_mb"]["rseqc"],
    conda:
        "envs/rseqc.yaml"
    log:
        f"{LOG}/rseqc/{{sample}}.bam_stat.log",
    shell:
        "bam_stat.py -i {input.bam} > {output} 2> {log}"

# b. Bibliotekos grandinės specifiškumas (strandedness)
rule rseqc_infer_experiment:
    input:
        bam=f"{OUT}/aligned/{{sample}}.sorted.bam",
        bed=GENE_BED,
    output:
        f"{OUT}/rseqc/{{sample}}.infer_experiment.txt",
    resources:
        mem_mb=config["memory_mb"]["rseqc"],
    conda:
        "envs/rseqc.yaml"
    log:
        f"{LOG}/rseqc/{{sample}}.infer_experiment.log",
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

# c. Skaitymų pasiskirstymas pagal genomo elementus
rule rseqc_read_distribution:
    input:
        bam=f"{OUT}/aligned/{{sample}}.sorted.bam",
        bed=GENE_BED,
    output:
        f"{OUT}/rseqc/{{sample}}.read_distribution.txt",
    resources:
        mem_mb=config["memory_mb"]["rseqc"],
    conda:
        "envs/rseqc.yaml"
    log:
        f"{LOG}/rseqc/{{sample}}.read_distribution.log",
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

# d (pagalbinė). Transkriptų poaibis gene body coverage skaičiavimui
rule gene_body_subset:
    input:
        GENE_BED,
    output:
        GENE_BODY_BED,
    params:
        step=config["rseqc"]["gene_body_subsample"],
    log:
        f"{LOG}/reference/gene_body_subset.log",
    shell:
        "awk -v n={params.step} 'NR % n == 0' {input} > {output} 2> {log}"

# d. Genų kūno padengimas (gene body coverage)
rule rseqc_gene_body_coverage:
    input:
        bam=f"{OUT}/aligned/{{sample}}.sorted.bam",
        bai=f"{OUT}/aligned/{{sample}}.sorted.bam.bai",
        bed=GENE_BODY_BED,
    output:
        txt=f"{OUT}/rseqc/{{sample}}.geneBodyCoverage.txt",
        pdf=f"{OUT}/rseqc/{{sample}}.geneBodyCoverage.curves.pdf",
    params:
        prefix=f"{OUT}/rseqc/{{sample}}",
    resources:
        mem_mb=config["memory_mb"]["rseqc"],
    conda:
        "envs/rseqc.yaml"
    log:
        f"{LOG}/rseqc/{{sample}}.gene_body_coverage.log",
    shell:
        "geneBody_coverage.py -r {input.bed} -i {input.bam} -o {params.prefix} &> {log}"

# e. Vidinis atstumas tarp porinių skaitymų (fragmento ilgio įvertinimui)
rule rseqc_inner_distance:
    input:
        bam=f"{OUT}/aligned/{{sample}}.sorted.bam",
        bed=GENE_BED,
    output:
        freq=f"{OUT}/rseqc/{{sample}}.inner_distance_freq.txt",
        pdf=f"{OUT}/rseqc/{{sample}}.inner_distance_plot.pdf",
    params:
        prefix=f"{OUT}/rseqc/{{sample}}",
    resources:
        mem_mb=config["memory_mb"]["rseqc"],
    conda:
        "envs/rseqc.yaml"
    log:
        f"{LOG}/rseqc/{{sample}}.inner_distance.log",
    shell:
        "inner_distance.py -i {input.bam} -r {input.bed} -o {params.prefix} &> {log}"

# f, g. Mėginių koreliacija ir PCA (deepTools)
rule multibamsummary:
    input:
        bams=BAMS,
        bais=BAIS,
    output:
        npz=f"{OUT}/deeptools/multiBamSummary.npz",
        counts=f"{OUT}/deeptools/multiBamSummary.counts.tab",
    params:
        labels=" ".join(samples),
        bs=config["deeptools"]["bin_size"],
    threads: config["threads"]["deeptools"]
    resources:
        mem_mb=config["memory_mb"]["deeptools"],
    conda:
        "envs/deeptools.yaml"
    log:
        f"{LOG}/deeptools/multiBamSummary.log",
    shell:
        """
        multiBamSummary bins -b {input.bams} --labels {params.labels} \
                        -bs {params.bs} -p {threads} \
                        -o {output.npz} --outRawCounts {output.counts} &> {log}
        """

rule plot_correlation:
    input:
        f"{OUT}/deeptools/multiBamSummary.npz",
    output:
        png=f"{OUT}/deeptools/correlation_heatmap.png",
        matrix=f"{OUT}/deeptools/correlation_matrix.tab",
    params:
        method=config["deeptools"]["cor_method"],
    conda:
        "envs/deeptools.yaml"
    log:
        f"{LOG}/deeptools/plotCorrelation.log",
    shell:
        """
        plotCorrelation -in {input} --corMethod {params.method} \
                        --whatToPlot heatmap --plotNumbers --skipZeros \
                        -o {output.png} --outFileCorMatrix {output.matrix} &> {log}
        """

rule plot_pca:
    input:
        f"{OUT}/deeptools/multiBamSummary.npz",
    output:
        png=f"{OUT}/deeptools/pca.png",
        data=f"{OUT}/deeptools/pca_data.tab",
    conda:
        "envs/deeptools.yaml"
    log:
        f"{LOG}/deeptools/plotPCA.log",
    shell:
        "plotPCA -in {input} -o {output.png} --outFileNameData {output.data} &> {log}"
