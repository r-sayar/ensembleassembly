rule fastqc_raw:
    input:
        fq=lambda wildcards: samples.at[wildcards.sample, "fq1" if wildcards.idx == 1 else "fq2"]
    output:
        html="results/qc/fastqc/{sample}_{idx}.html",
        zip="results/qc/fastqc/{sample}_{idx}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "results/logs/{sample}_{idx}_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v5.8.3/bio/fastqc"


rule fastqc_trim:
    input:
        fq=lambda wildcards: f"results/trim/{wildcards.sample}_{'forward_paired' if wildcards.idx == 1 else 'reverse_paired'}.fq.gz"
    output:
        html="results/qc/fastqc/{sample}_{idx}_trim.html",
        zip="results/qc/fastqc/{sample}_{idx}_fastqc_trim.zip"
    params:
        extra = "--quiet"
    log:
        "results/logs/{sample}_{idx}_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v5.8.3/bio/fastqc"