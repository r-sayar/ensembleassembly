checkpoint select_best_reference:
    input:
        multiqc_report="results/qc/multiqc_all/denovo/multiqc_data"#actually we need the /multiqc_data.json in this folder
    output:
        best_reference="results/best_reference.txt"
    script:
        "../scripts/select_best_reference.py"


def use_best_reference(wildcards):
    checkpoint_output = checkpoints.select_best_reference.get()
    with open(checkpoint_output.output.best_reference) as f:
        best_reference = "results/assembly/denovo/all_assemblies/"+f.readline().strip().split(": ")[1] + "_contigs.fasta"
    return best_reference

rule process_best_reference:
    input:
        best_reference=lambda wildcards: use_best_reference(wildcards)
    output:
        reference="results/references/{sample}_reference.fasta"
    shell:
        """
        cp {input.best_reference} {output.reference}
        """