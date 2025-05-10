#running this command isolated:
#snakemake --cores all --configfile config/config.yaml -s workflow/rules/read_classification.smk results/cleaned/sample1_classified_nonhuman_1.fq results/cleaned/sample1_classified_nonhuman_2.fq
#env needs to be activated before

rule classify_reads:
    # input:
    #     r1="/buffer/ag_bsc/pmsb_workflows_2025/team4_ensemble_assembly/DATA/SET1/{sample}_Illumina_MiSeq_paired_end_sequencing_1.fastq.gz",
    #     r2="/buffer/ag_bsc/pmsb_workflows_2025/team4_ensemble_assembly/DATA/SET1/{sample}_Illumina_MiSeq_paired_end_sequencing_2.fastq.gz"
    conda:
        "../env/krakentools.yaml"
    input:
        r1="results/trim/{sample}_trim_forward_paired.fq.gz",
        r2="results/trim/{sample}_trim_reverse_paired.fq.gz",     
    output:
        #classified and unclassified reads must have the same path! and need to set in shell as well
        #TODO: remove the Kraken2 from the name 
        kraken2_classified1="results/kraken2/{sample}/classified-reads_1.fq", 
        kraken2_unclassified1="results/kraken2/{sample}/unclassified-reads_1.fq",
        kraken2_classified2="results/kraken2/{sample}/classified-reads_2.fq",
        kraken2_unclassified2="results/kraken2/{sample}/unclassified-reads_2.fq",
        kraken2_report="results/kraken2/{sample}/{sample}.k2report",
        kraken2_stdout="results/kraken2/{sample}/{sample}.kraken2.stdout"
    params:
        kraken_db=config['kraken-db'],
    log:
        "results/logs/read_classification/{sample}-classify-reads.log",
    threads: 8
    shell:
        #kraken2 doesnt throw error if the input doesnt exist!! -> 
        """
        mkdir -p "results/kraken2/{wildcards.sample}"
        mkdir -p "results/logs/read_classification"

        kraken2 --db {params.kraken_db} \
            --gzip-compressed \
            --unclassified-out "results/kraken2/{wildcards.sample}/unclassified-reads#.fq" \
            --classified-out "results/kraken2/{wildcards.sample}/classified-reads#.fq" \
            --output {output.kraken2_stdout} \
            --report {output.kraken2_report} \
            --threads {threads} \
            --paired \
            {input.r1} {input.r2} 2> {log}
        """
        #might want to reduce / adjust threads usage 




rule clean_reads_with_kraken_test:
    conda:
        "../env/krakentools.yaml"
    input:
        # Input the paired-end classified files
        classified_r1="results/kraken2/{sample}/classified-reads_1.fq",
        classified_r2="results/kraken2/{sample}/classified-reads_2.fq",
        # *** ADD THIS: Input the standard Kraken2 output ***
        kraken2_stdout="results/kraken2/{sample}/{sample}.kraken2.stdout"
    output:
        # Output paired-end cleaned files
        cleaned_r1="results/cleaned/{sample}_classified_nonhuman_1.fq",
        cleaned_r2="results/cleaned/{sample}_classified_nonhuman_2.fq"
    log:
        "results/logs/read_classification/{sample}-clean-classified-reads-with-kraken-tools.log",
    params:
        human_taxid="9606"  # Taxonomic ID for humans #replace with the one from the config file to add more than one put space between them
    shell:
        """
        # Use KrakenTools to filter out human reads
        extract_kraken_reads.py \
            -k {input.kraken2_stdout} \
            -o {output.cleaned_r1} \
            -o2 {output.cleaned_r2} \
            -s1 {input.classified_r1} \
            -s2 {input.classified_r2} \
            --exclude --taxid {params.human_taxid} \
            --fastq-output \
            2> {log}
        """
