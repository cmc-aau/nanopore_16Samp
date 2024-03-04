rule qfilter:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}.fastq")
    output:
        temp(
            os.path.join(
                os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}_filtered.fastq")
            )
        ),
    log:
        os.path.join(config["log_dir"], "qfilter", "{sample}.log"),
    message:
        "Filtering reads using Filtlong"
    conda:
        "../envs/qfilter.yml"
    resources:
        mem_mb=512,
        runtime=10
    params:
        filtlong_args=config["filtlong_args"],
    threads: 1
    shell:
        """
        filtlong {params.filtlong_args} {input} > {output}
        """
