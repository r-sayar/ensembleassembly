INDEX_PREFIX = "results/index/human_index"

rule create_index:
    conda: "../env/mapping_env.yaml" 
    input:
        fasta=config["human_genome_fasta"]
    output:
        multiext(INDEX_PREFIX, ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "results/logs/create_index/bwa_index.log"
    threads: 8
    params:
        index_prefix=INDEX_PREFIX
    shell:
        """
        mkdir -p $(dirname {params.index_prefix})
        bwa index -p {params.index_prefix} {input.fasta} > {log} 2>&1       
        """


rule align_to_human_genome:
    conda: "../env/mapping_env.yaml" 
    input:
        cleaned_r1="results/cleaned/{sample}_classified_nonhuman_1.fq",
        cleaned_r2="results/cleaned/{sample}_classified_nonhuman_2.fq",
        index_files=multiext(INDEX_PREFIX, ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam=temp("results/aligned_host/{sample}.bam") 
    log:
        "results/logs/align_to_human_genome/{sample}.log"
    threads: 8 
    shell:
        """
        bwa mem -t {threads} {INDEX_PREFIX} {input.cleaned_r1} {input.cleaned_r2} | samtools view -bS - > {output.bam} 2> {log}
        """

rule remove_human_reads:
    conda:
        "../env/filtering_env.yaml" 
    input:
        bam="results/aligned_host/{sample}.bam"
    output:
        filtered_r1="results/filtered/{sample}_1.fq",
        filtered_r2="results/filtered/{sample}_2.fq",
    log:
        "results/logs/remove_human_reads/{sample}.log"
    shell:
        """
        samtools view -b -f 12 -F 256 {input.bam} | bedtools bamtofastq -i - -fq {output.filtered_r1} -fq2 {output.filtered_r2} 2> {log}
        """