#this file is unused as of now but might use later 
mkdir -p results/alignment results/consensus results/logs/consensus

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
