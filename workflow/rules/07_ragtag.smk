rule ragtag_fill_gaps:
    conda:
        "../env/ragtag_env.yaml"
        
    input:
        assembly="results/assembly/denovo/{sample}/contigs.fasta", # The assembly to be patched
        reference="results/references/{sample}_reference.fasta"  # The reference to use for patching
    output:
        # Directory where RagTag will write its intermediate patch outputs
        ragtag_patch_dir=temp(directory("results/ragtag_workdir/{sample}_patched")),
        # The final patched assembly file we want to keep
        patched_assembly="results/assembly/patched_reference/{sample}_patched_contigs.fasta" # Changed output path for clarity
    log:
        "logs/ragtag_patch_{sample}.log" # Adjusted log file name
    shell:
        """
        # Run ragtag patch, specifying the output directory
        ragtag.py patch {input.reference} {input.assembly} -o {output.ragtag_patch_dir} > {log} 2>&1
        
        # RagTag patch typically creates 'ragtag.patch.fasta' in its output directory.
        # Move this specific file to the desired final output path and name.
        mv {output.ragtag_patch_dir}/ragtag.patch.fasta {output.patched_assembly}
        """