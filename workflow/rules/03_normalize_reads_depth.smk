rule normalize_reads_depth:
    conda:
        "../env/bbmap_env.yaml"
    input:
        r1="results/filtered/{sample}_1.fq",
        r2="results/filtered/{sample}_2.fq",
    output:
        norm_r1="results/normalized/{sample}_normalized_1.fq",
        norm_r2="results/normalized/{sample}_normalized_2.fq",
        hist="results/normalized/{sample}_depth_histogram.txt"
    params:
        target_depth=config["normalize_reads"]["target_depth"],
        min_depth=config["normalize_reads"]["min_depth"]
    log:
        "results/logs/normalize_reads_depth/{sample}_normalize.log"
    threads: 8
    shell:
        """
        bbnorm.sh \
            in1={input.r1} \
            in2={input.r2} \
            out1={output.norm_r1} \
            out2={output.norm_r2} \
            target={params.target_depth} \
            min={params.min_depth} \
            hist={output.hist} \
            threads={threads} \
            overwrite=t \
            2> {log}
        """