rule generate_consensus:
    input:
        reference_assembly="results/{sample}_{reference}_filled.fasta",
        r1="results/trim/{sample}_trim_forward_paired.fq.gz",
        r2="results/trim/{sample}_trim_reverse_paired.fq.gz"
    output:
        consensus_fasta="results/consensus/{sample}_{reference}_consensus.fasta",
    params:
        alignment_bam="results/alignment/{sample}_{reference}.bam", # Intermediate BAM
        sorted_bam="results/alignment/{sample}_{reference}.sorted.bam", # Path for sorted BAM
        vcf_file="results/consensus/{sample}_{reference}.vcf.gz" # Path for VCF
    log:
        "results/logs/consensus/{sample}_{reference}_consensus.log"
    conda:
        "../env/bwa_samtools_bcftools.yaml" 
    threads: 8
    shell:
        #do we need to add the "final" from stages here?
        """
        mkdir -p results/alignment results/consensus results/logs/consensus

        echo "Starting consensus generation for {wildcards.sample} against {wildcards.reference}" > {log}

        bwa index {input.reference_assembly} >> {log} 2>&1

        bwa mem -t {threads} {input.reference_assembly} {input.r1} {input.r2} 2>> {log} | \
            samtools view -Sb -F 4 -@ {threads} -o {params.alignment_bam} - >> {log} 2>&1

        samtools sort -@ {threads} -o {params.sorted_bam} {params.alignment_bam} >> {log} 2>&1

        samtools index -@ {threads} {params.sorted_bam} >> {log} 2>&1

        samtools mpileup -uf {input.reference_assembly} {params.sorted_bam} 2>> {log} | \
            bcftools call -c --ploidy 1 -Oz -o {params.vcf_file} >> {log} 2>&1

        bcftools index {params.vcf_file} >> {log} 2>&1

        bcftools consensus -f {input.reference_assembly} {params.vcf_file} -o {output.consensus_fasta} >> {log} 2>&1

        echo "Consensus generation for {wildcards.sample} against {wildcards.reference} finished." >> {log}

        """

#output has to be in the format:         
    #"results/assembly/denovo/all_assemblies/{sample}_contigs.fasta"
