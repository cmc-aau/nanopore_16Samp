rule qfilter:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}.fastq")
    output:
        fastq = temp(
            os.path.join(
                os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}_filtered.fastq")
            )
        ),
        totalfilteredreads_file = temp(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "{sample}_totalfilteredreads.csv"
            )
        ),
    log:
        os.path.join(config["log_dir"], "qfilter", "{sample}.log"),
    message:
        "Filtering reads using Filtlong"
    conda:
        "../envs/qfilter.yml"
    resources:
        mem_mb=lambda wc, input: max(3 * input.size_mb, 512),
        runtime=30
    params:
        filtlong_args=config["filtlong_args"],
    threads: 1
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail
        
        filtlong {params.filtlong_args} {input} > {output.fastq}

        # calc total number of reads
        num_reads=$(($(wc -l < "{output.fastq}") / 4))
        echo "{wildcards.sample},$num_reads" > "{output.totalfilteredreads_file}"
        """
