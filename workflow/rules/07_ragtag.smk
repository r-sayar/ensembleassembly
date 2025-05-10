rule ragtag_fill_gaps:
    conda:
        "envs/ragtag.yaml"
    input:
        assembly="results/references/{sample}_{reference}.fasta",
        reference="results/{reference}.fasta"
    output:
        filled_assembly="results/{sample}_{reference}_filled.fasta"
    log:
        "logs/ragtag_fill_gaps_{sample}_{reference}.log"
    shell:
        """
        ragtag.py fillgap {input.reference} {input.assembly} -o results/{wildcards.sample}_{wildcards.reference}_filled > {log} 2>&1
        mv results/{wildcards.sample}_{wildcards.reference}_filled/ragtag.filled.fasta {output.filled_assembly}
        """