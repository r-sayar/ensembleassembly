rule fastqc_raw:
    input:
        fq=lambda wildcards: samples.at[wildcards.sample, "fq1" if wildcards.idx == '1' else "fq2"]
    output:
        html="results/qc/fastqc/raw/{sample}_{idx}.html",
        zip="results/qc/fastqc/raw/{sample}_{idx}_fastqc.zip"
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
        html="results/qc/fastqc/trimmed/{sample}_trim_{idx}.html",
        zip="results/qc/fastqc/trimmed/{sample}_trim_{idx}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "results/logs/{sample}_{idx}_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v5.8.3/bio/fastqc"


rule multiqc_raw:
    input:
        expand("results/qc/fastqc/raw/{sample}_{idx}.html", idx=['1','2'], sample=samples.index)  # Raw reads only
    output:
        report="results/qc/multiqc_raw/multiqc_report.html",
        report_data=directory("results/qc/multiqc_raw/multiqc_data/")
    conda:
        "../env/myenv.yaml"
    log:
        "results/logs/multiqc_raw.log"
    shell:
        """
        multiqc results/qc/fastqc/raw -o results/qc/multiqc_raw > {log} 2>&1
        """


rule multiqc_trimmed:
    input:
        expand("results/qc/fastqc/trimmed/{sample}_trim_{idx}.html", idx=['1','2'], sample=samples.index)  # Trimmed reads only
    output:
        report="results/qc/multiqc_trimmed/multiqc_report.html",
        report_data=directory("results/qc/multiqc_trimmed/multiqc_data/")
    conda:
        "../env/myenv.yaml"
    log:
        "results/logs/multiqc_trimmed.log"
    shell:
        """
        multiqc results/qc/fastqc/trimmed -o results/qc/multiqc_trimmed > {log} 2>&1
        """


rule multiqc_all:
    input:
        # FastQC
        expand("results/qc/fastqc/raw/{sample}_{idx}.html", idx=['1','2'], sample=samples.index),
        expand("results/qc/fastqc/trimmed/{sample}_trim_{idx}.html", idx=['1','2'], sample=samples.index),
        # BUSCO + QUAST (denovo and final)
        lambda wildcards: expand(f"results/qc/quast/{wildcards.stage}/{{sample}}/report.tsv", sample=samples.index),
        lambda wildcards: expand(f"results/qc/busco/{wildcards.stage}/reports/short_summary_{wildcards.stage}_{{sample}}.txt", sample=samples.index),
        # Kraken2
        expand("results/kraken2/{sample}/{sample}.k2report", sample=samples.index)

    output:
        report="results/qc/multiqc_all/{stage}/multiqc_report.html",
        report_data=directory("results/qc/multiqc_all/{stage}/multiqc_data/")
        
    conda:
        "../env/myenv.yaml"
    log:
        "results/logs/multiqc-{stage}.log"
    shell:
        """
        #maybe this is an issue 
        #mkdir -p results/qc/quast/{wildcards.stage}
        #mkdir -p results/qc/busco/{wildcards.stage}/reports
        #mkdir -p results/qc/multiqc_all/{wildcards.stage}

        multiqc results/qc/fastqc results/qc/quast/{wildcards.stage} results/qc/busco/{wildcards.stage}/reports results/kraken2 -o results/qc/multiqc_all/{wildcards.stage} > {log} 2>&1
        """


# (separat aufrufbar)
rule busco_download:
    output:
        lineage_dir = directory("busco_downloads/lineages/{lineage}")
    conda:
        "../env/busco.yaml"
    shell:
        """
        busco --download {wildcards.lineage} -o temp_download -l {output.lineage_dir}
        """


rule busco:
    input:
        assembly="results/assembly/{stage}/all_assemblies/{sample}_contigs.fasta",
        lineage_dir = directory(f"busco_downloads/lineages/{config['busco_lineage']}")
    output:
        f"results/qc/busco/{{stage}}/{{sample}}/short_summary.specific.{config['busco_lineage']}.{{sample}}.txt"
    conda:
        "../env/busco.yaml"
    threads: 5
    log:
        "results/logs/busco/{stage}/{sample}.log"
    shell:
        """
        busco -q -c {threads} -f -m genome -l {input.lineage_dir} -o results/qc/busco/{wildcards.stage}/{wildcards.sample} -i {input.assembly} > {log} 2>&1
        """


rule busco_summary_for_multiqc:
    input:
        f"results/qc/busco/{{stage}}/{{sample}}/short_summary.specific.{config['busco_lineage']}.{{sample}}.txt"
    output:
        "results/qc/busco/{stage}/reports/short_summary_{stage}_{sample}.txt"
    shell:
        """
        cp {input} {output}
        """


rule quast:
    input:
        assembly="results/assembly/{stage}/all_assemblies/{sample}_contigs.fasta"
    output:
        report = "results/qc/quast/{stage}/{sample}/report.tsv"
    conda:
        "../env/quast.yaml"
    threads: 5
    log:
        "results/logs/quast/{stage}/{sample}/quast.log",
    shell:
        """
        quast -t {threads} --glimmer -o results/qc/quast/{wildcards.stage}/{wildcards.sample} --labels {wildcards.stage}_{wildcards.sample} {input.assembly} > {log} 2>&1
        """