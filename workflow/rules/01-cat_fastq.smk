import glob
import os


rule concatenate_fastq:
    input:
        # function to list all fastq files per wildcard (subfolder/sample)
        # see https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#input-functions
        lambda wildcards: glob.glob(
            os.path.join(config["input_dir"], wildcards.sample, "*.fastq.gz")
        ),
    output:
        fastq=temp(
            os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}.fastq")
        ),
        total_reads_file=temp(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "{sample}_totalreads.csv"
            )
        ),
    log:
        os.path.join(config["log_dir"], "concatenate_fastq", "{sample}.log"),
    resources:
        mem_mb=512,
        runtime=30
    conda:
        "../envs/gzip.yml"
    threads: 1
    message:
        "Decompressing reads if gzip'ed and concatenating into a single file, then calculate total number of reads"
    shell:
        """
        # decompress only if compressed, but concatenate regardless
        gunzip -cdfq {input} > {output.fastq}
        
        # calc total number of reads
        num_reads=$(($(wc -l < "{output.fastq}") / 4))
        echo "{wildcards.sample},$num_reads" > "{output.total_reads_file}"
        """
