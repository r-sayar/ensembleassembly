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


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_{idx}.html", idx=['1','2'], sample=samples.index),
        expand("results/qc/fastqc/{sample}_{idx}_trim.html", idx=['1','2'], sample=samples.index)
    output:
        report="results/qc/multiqc_report.html",
        report_data = directory("results/qc/multiqc_data/")
    conda:
        "../env/myenv.yaml"
    log:
        "results/logs/multiqc.log"
    shell:
        """
        multiqc results/qc/fastqc/ -o results/qc > {log} 2>&1
        """