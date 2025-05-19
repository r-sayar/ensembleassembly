rule ragtag_fill_gaps:
    conda:
        "envs/ragtag.yaml"
    input:
        assembly="results/assembly/denovo/{sample}/contigs.fasta",
        reference="results/references/{sample}_reference.fasta"
    output:
        ragtag=temp("results/{sample}_reference_filled"),
        filled_assembly="results/{sample}_reference_filled.fasta"
    log:
        "logs/ragtag_fill_gaps_{sample}_reference.log"
    shell:
        """
        ragtag.py fillgap {input.reference} {input.assembly} -o {output.ragtag} > {log} 2>&1
        mv {output.ragtag} {output.filled_assembly}
        """