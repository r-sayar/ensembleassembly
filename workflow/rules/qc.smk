rule fastqc_raw:
    input:
        fq=lambda wildcards: samples.at[wildcards.sample, "fq1" if wildcards.idx == '1' else "fq2"]
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
        fq=lambda wildcards: f"results/trim/{wildcards.sample}_trim_{'forward_paired' if wildcards.idx == '1' else 'reverse_paired'}.fq.gz"
    output:
        html="results/qc/fastqc/{sample}_trim_{idx}.html",
        zip="results/qc/fastqc/{sample}_trim_{idx}_fastqc.zip"
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
        expand("results/qc/fastqc/{sample}_trim_{idx}.html", idx=['1','2'], sample=samples.index),
        expand("results/qc/quast/{sample}/report.html", sample=samples.index),
        expand("results/qc/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.tsv", sample=samples.index),
        "results/qc/busco/all_stats.tsv"
    output:
        report="results/qc/multiqc_report.html",
        report_data = directory("results/qc/multiqc_data/")
    conda:
        "../env/myenv.yaml"
    log:
        "results/logs/multiqc.log"
    shell:
        """
        multiqc results/qc -c config/my_multiqc_config.yaml -o results/qc > {log} 2>&1
        """


rule busco:
    input:
        assembly="results/assembly/{sample}/contigs.fasta"
    output:
        "results/qc/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.txt"
    conda:
        "../env/busco.yaml"
    threads: 5
    log:
        "results/logs/busco/{sample}.log"
    shell:
        """
        busco -q -c {threads} -f -m genome -l bacteria -o results/qc/busco/{wildcards.sample} -i {input.assembly} > {log} 2>&1
        """


rule busco2tsv:
    threads: 1
    input:
        "results/qc/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.txt"
    output:
        "results/qc/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.tsv"
    shell:
        r"""
        perl -lne 'print "{wildcards.sample}\t$1\t$2\t$3\t$4" if /C:([\-\d\.]+).*F:([\-\d\.]+).*M:([\-\d\.]+).*n:(\d+)/' {input} > {output}
        """


rule gather_stats_busco:
    input:
        expand("results/qc/busco/{sample}/short_summary.specific.bacteria_odb10.{sample}.tsv", sample=samples.index)
    output:
        "results/qc/busco/all_stats.tsv"
    shell:
        """
        cat {input} > {output}
        """




rule quast:
    input:
        assembly="results/assembly/{sample}/contigs.fasta"
    output:
        "results/qc/quast/{sample}/report.html"
    conda:
        "../env/quast.yaml"
    threads: 5
    log:
        "results/logs/quast/{sample}/quast.log",
    shell:
        """
        quast -t {threads} --glimmer -o results/qc/quast/{wildcards.sample} {input.assembly} > {log} 2>&1
        """

