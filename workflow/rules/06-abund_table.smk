rule abund_table:
    input:
        totalreads_files=expand(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "{sample}_totalreads.csv"
            ),
            sample=sample_dirs,
        ),
        totalfilteredreads_files=expand(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "{sample}_totalfilteredreads.csv"
            ),
            sample=sample_dirs,
        ),
        mapping_overviews=expand(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "{sample}.idmapped.txt"
            ),
            sample=sample_dirs,
        ),
    output:
        os.path.join(config["output_dir"], "summary.txt"),
        os.path.join(config["output_dir"], "otutable_mappedreads.tsv"),
        os.path.join(config["output_dir"], "otutable_normalised.tsv"),
    log:
        os.path.join(config["log_dir"], "abund_table.log"),
    resources:
        mem_mb=32768,
        runtime=60
    message:
        "Generating abundance table and writing final output files"
    threads: min(config["max_threads"], 4)
    conda:
        "../envs/abund_table.yml"
    script:
        "../scripts/abund_table.R"
