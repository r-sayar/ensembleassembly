rule spades_assembly:
    input:
        fq1 = lambda wildcards: f"results/trim/{wildcards.sample}_trim_forward_paired.fq.gz",
        fq2 = lambda wildcards: f"results/trim/{wildcards.sample}_trim_reverse_paired.fq.gz"
    output:
        contigs = "results/assembly/{sample}/contigs.fasta",
        scaffolds = "results/assembly/{sample}/scaffolds.fasta"
    conda:
        "../env/myenv.yaml"
    log:
        "results/logs/{sample}_assembly.log"
    threads: 8
    shell:
        """
        spades.py \
            -1 {input.fq1} \
            -2 {input.fq2} \
            -o results/assembly/{wildcards.sample} \
            -t {threads} \
            --only-assembler \
            > {log} 2>&1
        """


rule copy_assemblies:
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/assembly/all_assemblies/{sample}_contigs.fasta"
    shell:
        """
        cp {input} {output}
        """
