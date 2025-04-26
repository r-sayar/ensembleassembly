rule normalize_reads_depth:
    input:
        # Adjust path based on previous rule output
        r1="results/cleaned/{sample}_classified_nonhuman_1.fq",
        r2="results/cleaned/{sample}_classified_nonhuman_2.fq"
    output:
        # Define output paths for normalized reads
        norm_r1="results/normalized/{sample}_normalized_1.fq",
        norm_r2="results/normalized/{sample}_normalized_2.fq",
        hist="results/normalized/{sample}_depth_histogram.txt" # Optional histogram output
    params:
        target_depth=100,  # Set the desired target depth (adjust as needed)
        min_depth=5      # Optional: minimum depth to keep reads
    log:
        "results/logs/normalize_reads_depth/{sample}_normalize.log"
    conda:
        "../env/bbmap_env.yaml"  # Path to your conda environment YAML file for bbmap
    threads: 8 # bbnorm can use multiple threads
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