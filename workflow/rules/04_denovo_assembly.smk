rule spades_assembly:
    input:
        fq1="results/normalized/{sample}_normalized_1.fq",
        fq2="results/normalized/{sample}_normalized_2.fq",
    output:
        contigs = "results/assembly/denovo/{sample}/contigs.fasta",
        scaffolds = "results/assembly/denovo/{sample}/scaffolds.fasta"
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
            -o results/assembly/denovo/{wildcards.sample} \
            -t {threads} \
            --only-assembler \
            > {log} 2>&1
        """


rule copy_assemblies:
    input:
        "results/assembly/denovo/{sample}/contigs.fasta"
    output:
        "results/assembly/denovo/all_assemblies/{sample}_contigs.fasta"
    shell:
        """
        cp {input} {output}
        """
