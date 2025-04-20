rule classify_reads:
    input:
        r1="/buffer/ag_bsc/pmsb_workflows_2025/team4_ensemble_assembly/DATA/SET1/{sample}_Illumina_MiSeq_paired_end_sequencing_1.fastq.gz",
        r2="/buffer/ag_bsc/pmsb_workflows_2025/team4_ensemble_assembly/DATA/SET1/{sample}_Illumina_MiSeq_paired_end_sequencing_2.fastq.gz"
    output:
        #classified and unclassified reads must have the same path! and need to set in shell as well
        #TODO: remove the Kraken2 from the name 
        kraken2_classified1="results/kraken2/{sample}/classified-reads_1.fq.gz", 
        kraken2_unclassified1="results/kraken2/{sample}/unclassified-reads_1.fq.gz",
        kraken2_classified2="results/kraken2/{sample}/classified-reads_2.fq.gz",
        kraken2_unclassified2="results/kraken2/{sample}/unclassified-reads_2.fq.gz",
        kraken2_report="results/kraken2/{sample}/{sample}.k2report",
        kraken2_stdout="results/kraken2/{sample}/{sample}.kraken2.stdout"
    params:
        kraken_db=config['kraken-db'],
    log:
        "results/logs/read_classification/{sample}-classify-reads.log",
    threads: 8
    shell:
        """
        kraken2 --db {params.kraken_db} \
            --unclassified-out "results/kraken2/{wildcards.sample}/unclassified-reads#.fq.gz" \ 
            --classified-out "results/kraken2/{wildcards.sample}/classified-reads#.fq.gz" \
            --output {output.kraken2_stdout} \
            --report {output.kraken2_report} \
            --threads {threads} \
            --paired \
            {input.r1} {input.r2} 2> {log}
        """
        #might want to reduce / adjust threads usage 
rule clean_reads:
    conda:
        "../env/bbmap_env.yaml"
    input:
        # Input the paired-end classified files
        classified_r1="results/kraken2/{sample}/classified-reads_1.fq.gz",
        classified_r2="results/kraken2/{sample}/classified-reads_2.fq.gz",
        # Input the standard Kraken2 output (needed for awk)
        kraken2_stdout="results/kraken2/{sample}/{sample}.kraken2.stdout"
    output:
        # Output paired-end cleaned files
        cleaned_r1="results/cleaned/{sample}_classified_nonhuman_1.fq.gz",
        cleaned_r2="results/cleaned/{sample}_classified_nonhuman_2.fq.gz"
    log:
        "results/logs/read_classification/{sample}-clean-classified-reads.log",
    params:
        human_taxid="9606"  # Taxonomic ID for humans
    shell:
        """
        # Temporary file for human read IDs
        HUMAN_IDS_TMP=$(mktemp --tmpdir={wildcards.sample}_human_read_ids.XXXXXX.txt)

        # Extract read IDs classified directly as human from the standard output
        # Kraken2 stdout format: C/U <tab> ReadID <tab> TaxID ...
        # Extract the base ReadID (removing potential /1 or /2 suffix)
        awk -v taxid={params.human_taxid} \
            'BEGIN {{FS="\\t"; OFS="\\t"}} $3 == taxid {{ split($2, id_part, "/"); print id_part[1] }}' \
            {input.kraken2_stdout} > "$HUMAN_IDS_TMP"

        # Filter the *paired-end classified* reads using BBTools filterbyname.sh
        # Use 'in1=...', 'in2=...' for paired input, 'out1=...', 'out2=...' for paired output
        # include=f removes reads listed in the names file
        # Assumes filterbyname.sh is in PATH
        filterbyname.sh \
            in1={input.classified_r1} \
            in2={input.classified_r2} \
            out1={output.cleaned_r1} \
            out2={output.cleaned_r2} \
            names="$HUMAN_IDS_TMP" \
            include=f \
            overwrite=t \
            qin=auto qout=auto \
            2> {log} # Capture stderr (filterbyname logs here)

        # Clean up temporary file
        rm "$HUMAN_IDS_TMP"
        """

#running this command isolated:
#snakemake --cores all --configfile config/config.yaml -s workflow/rules/read_classification.smk results/cleaned/sample1_classified_nonhuman_interleaved.fq.gz 
#env needs to be activated before