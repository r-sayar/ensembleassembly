rule trimm:
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        forward_paired="results/trim/{sample}_forward_paired.fq.gz",
        forward_unpaired="results/trim/{sample}_forward_unpaired.fq.gz",
        reverse_paired="results/trim/{sample}_reverse_paired.fq.gz",
        reverse_unpaired="results/trim/{sample}_reverse_unpaired.fq.gz"
    conda:
        "../env/myenv.yaml"
    log:
        "results/logs/{sample}_trimmomatic.log"
    shell:
        """
        trimmomatic PE \
        {input.r1} {input.r2} \
        {output.forward_paired} {output.forward_unpaired} \
        {output.reverse_paired} {output.reverse_unpaired} \
        ILLUMINACLIP:{config[trimmomatic][adapters]}:{config[trimmomatic][illuminaclip]} \
        LEADING:{config[trimmomatic][leading]} \
        TRAILING:{config[trimmomatic][trailing]} \
        SLIDINGWINDOW:{config[trimmomatic][slidingwindow]} \
        MINLEN:{config[trimmomatic][minlen]} \
        &> {log}
        """