# could have used the minimap2 and samtools wrappers,
# but piping directly to samtools view from minimap2 is most efficient
rule map2db:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}_filtered_renamed.fastq"),
    output:
        temp(os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}.sam")),
    log:
        os.path.join(config["log_dir"], "map2db", "{sample}.log"),
    resources:
        #depending on the tool memory usage usually scales with threads 
        #and/or input/database file size. Can calculate dynamically
        mem_mb=24000,
        runtime=60
    threads: config["max_threads"]
    message:
        "Mapping against database and filtering output"
    params:
        db_fasta=config["db_fasta"],
    conda:
        "../envs/map2db.yml"
    shell:
        """
    minimap2 \
      -ax map-ont \
      -K20M \
      -t {threads} \
      --secondary=no \
      {params.db_fasta} \
      {input} \
      | samtools \
        view \
        -F 4 \
        -F 256 \
        -F 2048 \
        --threads {threads} \
        -o {output}
    """
