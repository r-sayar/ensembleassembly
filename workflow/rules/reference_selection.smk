checkpoint select_best_reference:
    input:
        multiqc_report="results/qc/multiqc_all/multiqc_data/multiqc_data.json"
    output:
        best_reference="results/best_reference.txt"
    run:
        import json
        
        with open(input.multiqc_report) as f:
            multiqc_data = json.load(f)
        
        references = []
        for sample in multiqc_data.get("report_data_sources", {}).get("busco", []):
            busco_score = sample.get("Complete_BUSCOs (%)", 0) #issue potential if name changes
            quast_data = next((q for q in multiqc_data.get("report_data_sources", {}).get("quast", []) if q["Sample"] == sample["Sample"]), None)
            if quast_data:
                n50 = quast_data.get("N50", 0)
                references.append((sample["Sample"], busco_score, n50))
        
        best_reference = max(references, key=lambda x: (x[1], x[2]))
        
        with open(output.best_reference, "w") as f:
            f.write(f"Best Reference: {best_reference[0]}\n")
            f.write(f"BUSCO Score: {best_reference[1]}\n")
            f.write(f"N50: {best_reference[2]}\n")


def use_best_reference(wildcards):
    checkpoint_output = checkpoints.select_best_reference.get()
    with open(checkpoint_output.output.best_reference) as f:
        best_reference = f.readline().strip().split(": ")[1]
    return best_reference

rule process_best_reference:
    input:
        best_reference=lambda wildcards: use_best_reference(wildcards)
    output:
        references="results/references/{sample}_{reference}.fasta"
    shell:
        """
        cp data/references/{input.best_reference}.fasta {output.references}
        """