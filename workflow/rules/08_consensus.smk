rule generate_consensus:
    input:
        reference_assembly="results/assembly/patched_reference/{sample}_{reference}_patched_contigs.fasta" ,

        r1="results/trim/{sample}_trim_forward_paired.fq.gz",
        r2="results/trim/{sample}_trim_reverse_paired.fq.gz"
    output:
        #this needs to be something else
        assembly="results/assembly/final/all_assemblies/{sample}_{reference}_contigs.fasta"

    params:
        alignment_bam="results/alignment/{sample}_reference.bam", # Intermediate BAM
        sorted_bam="results/alignment/{sample}_reference.sorted.bam", # Path for sorted BAM
        vcf_file="results/consensus/{sample}_reference.vcf.gz" # Path for VCF
    log:
        "results/logs/consensus/{sample}_{reference}_consensus.log"
    conda:
        "../env/bwa_samtools_bcftools.yaml" 
    threads: 8
    shell:
        #do we need to add the "final" from stages here?
        #we also should split this up into several rules -> there can be a lot of intermediate files
        """
        mkdir -p results/alignment results/consensus results/logs/consensus

        bwa index {input.reference_assembly} >> {log} 2>&1

        bwa mem -t {threads} {input.reference_assembly} {input.r1} {input.r2} 2>> {log} | \
            samtools view -Sb -F 4 -@ {threads} -o {params.alignment_bam} - >> {log} 2>&1

        samtools sort -@ {threads} -o {params.sorted_bam} {params.alignment_bam} >> {log} 2>&1

        samtools index -@ {threads} {params.sorted_bam} >> {log} 2>&1

        samtools consensus {params.sorted_bam} -o {output.assembly}

        """

#output has to be in the format:         
    #"results/assembly/denovo/all_assemblies/{sample}_contigs.fasta"
