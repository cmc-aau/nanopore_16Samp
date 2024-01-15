import glob
import os


# function to list all fastq files per wildcard (subfolder/sample)
def listFastq(wildcards):
    fastqs = glob.glob(
        os.path.join(config["input_dir"], wildcards.sample, "*.fastq.gz")
    )
    return fastqs


rule concatenate_fastq:
    input:
        listFastq,
    output:
        fastq=temp(
            os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}.fastq")
        ),
        total_reads_file=os.path.join(
            config["tmp_dir"], "samples", "{sample}", "{sample}_totalreads.csv"
        ),
    resources:
        mem_mb=600,
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
