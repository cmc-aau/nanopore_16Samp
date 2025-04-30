rule mapping_overview:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}.sam"),
    output:
        temp(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "{sample}.idmapped.txt"
            )
        ),
    log:
        os.path.join(config["log_dir"], "mapping_overview", "{sample}.log"),
    message:
        "Creating mapping overview"
    conda:
        "../envs/checksum_fq_headers.yml"
    resources:
        mem_mb=1024,
        runtime=30
    threads: 1
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail
        
        sed ' /^@/ d' "{input}" | \
        awk '{{
            for(i=1;i<=NF;i++){{
            if($i ~ /^NM:i:/){{sub("NM:i:", "", $i); mm = $i}}
            }}
            split($6, count, /[^0-9]+/);
            split($6, type, /[^A-Z]*/);
            for(i=1; i<= length(count)-1; i++){{
            if(type[i + 1] ~ /[DIM]/){{aln+=count[i]}};
            }}
            print $1, $2, $3, length($10), aln, (aln - mm)/aln, $12, $14, $20
            aln=0;
        }}' \
        > "{output}"
        """
