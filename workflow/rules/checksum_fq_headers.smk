rule checksum_fq_headers:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}.fastq"),
    output:
        renamed_fastq=temp(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "{sample}_renamed.fastq"
            )
        ),
        seqid_checksums=temp(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}", "seqid_checksums.tsv"
            )
        ),
    log:
        os.path.join(config["log_dir"], "checksum_fq_headers", "{sample}.log"),
    conda:
        "../envs/checksum_fq_headers.yml"
    message:
        "Replacing sequence headers in fastq file with sha-256 checksums of themselves"
    resources:
        threads=1,
        mem_mb=lambda wc, input: max(0.5 * input.size_mb, 512),
        runtime=10
    threads: 1
    shell:
        """
      awk 'NR%4 == 1 {{print $0}}' "{input}" |\
        perl -MDigest::SHA=sha256_hex -nlE'say"$_\t".sha256_hex($_)' >\
        "{output.seqid_checksums}"

      awk 'FNR == NR {{
        header[FNR]="@"$NF
        next
        }}
        {{
          if (FNR%4 == 1) {{
            print header[(FNR-1)/4+1]
          }}
          if (FNR%4 != 1) {{
            print $0
          }}
        }}' \
        "{output.seqid_checksums}" \
        "{input}" >\
        "{output.renamed_fastq}"
    """
