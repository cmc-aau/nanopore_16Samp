rule abund_table:
    input:
        totalreads_files=expand(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "{sample}_totalreads.csv"
            ),
            sample=sample_dirs,
        ),
        mapping_overviews=expand(
            os.path.join(
                config["output_dir"], "samples", "{sample}", "{sample}.idmapped.txt"
            ),
            sample=sample_dirs,
        ),
    output:
        os.path.join(config["output_dir"], "mappings_detailed.txt"),
        os.path.join(config["output_dir"], "summary.txt"),
        os.path.join(config["output_dir"], "otutable_mappedreads.tsv"),
        os.path.join(config["output_dir"], "otutable_normalised.tsv"),
        os.path.join(config["output_dir"], "totalreads.csv"),
    resources:
        mem_mb=10240,
    message:
        "Generating abundance table and writing final output files"
    threads: min(config["max_threads"], 8)
    conda:
        "../envs/abund_table.yml"
    script:
        "../scripts/abund_table.R"
