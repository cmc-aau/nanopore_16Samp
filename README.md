# Snakemake workflow: `cmc-aau/nanopore_16Samp`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.18.2-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/cmc-aau/nanopore_16Samp/workflows/Tests/badge.svg)](https://github.com/cmc-aau/nanopore_16Samp/actions?query=branch%3Amain+workflow%3ATests)

"KISS" workflow to map any amplicon reads from the ONT platforms against a taxonomic database and produce an abundance table compatible with the [ampvis2](https://github.com/kasperskytte/ampvis2) R package. 

Requires a good database with high coverage of the diversity in the environment sampled for reliable results since no de-novo OTUs/ASVs are inferred from the reads (the nanopore platform still has too high error-rate for that), instead the sequences in the database itself are used and abundance estimated by counting mapping hits per sample/barcode. Expect some over-classification since we are not using classification algorithms but just mapping reads to a database. For activated sludge or anaerobic digester samples we have obtained really good results using the MiDAS ecosystem specific database (https://midasfieldguide.org). Comparing to Illumina data from the same samples we obtain almost identical results.

## Usage
The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=cmc-aau/nanopore_16Samp).

First install `snakemake` and `snakedeploy`, then deploy the workflow in your project using fx

```
snakedeploy deploy-workflow https://github.com/cmc-aau/nanopore_16Samp . --tag v123
```

Replace `v123` with one of the versions listed under releases.

Then adjust the configuration file `config.yml` (described in [`config/README.md`](config/README.md)) and run the workflow with `snakemake` according to your setup.

## Workflow rules/steps
Steps 1-5 are performed for each barcode individually and run in parallel (for each subfolder in the input folder)
1. `concatenate_fastq`: Decompress reads if gzip'ed and concatenate all files into a single file, then calculate total number of reads.
2. `qfilter`: Runs `filtlong` (usually just for length and Q-score filtering), then calculate total number of reads.
3. `checksum_fq_headers`: Replace sequence headers in the fastq files with sha-256 checksums of themselves (sometimes sequence headers are too long for some tools).
4. `map2db`: Map reads to database using `minimap2`, pipe through `samtools`, and output a SAM file.
5. `mapping_overview`: Count number of reads that mapped to a database sequence and calculate some stats to use for post-filtering.
6. `abund_table`: Filter mappings based on stats, and then boil it all down to a single abundance table per sample/barcode including taxonomy of the database reference.

## Output
3 files are output in the `config['output_dir']` directory:
 - `otutable_mappedreads.tsv`: Database matches, summed up per sample. Taxonomy of the database reference is appended in the last columns.
 - `otutable_normalised.tsv`: Database matches, summed up per sample, normalised to the total reads per sample BEFORE any filtering. Taxonomy of the database reference is appended in the last columns. You should use this one in most cases.
 - `summmary.txt`: Filtering and mapping stats per sample/barcode. Details below.

### `summary.txt`
Description of columns in the `summary.txt` file:

| Column | Description |
| --- | --- |
| `barcode` | Sample/barcode name. |
| `median_alnlen` | The median alignment length with the database reference of all mapping hits in the database. |
| `median_querylen` | The median length of all query sequences. |
| `median_IDfrac` | The median identity fraction between all mapping hits and the database sequences (`0.0-1.0`). |
| `mapped_reads` | Number of reads that mapped to a database sequence (= mapping hits) before any post-filtering of the mapping results. |
| `total_reads` | Total reads in the raw data before any pre-filtering. |
| `total_filtered_reads` | Total reads after pre-filtering (before mapping). |
| `pct_mapped` | Percent of pre-filtered reads that mapped to a database sequence (`mapped_reads / total_filtered_reads * 100`). |
| `nfilt` | Number of reads (mapping hits) after post-filtering (based on the `minalignlen` and `minIDfrac` values set in the config). |
| `pct_filt` | Percent of reads (mapping hits) that have been filtered (`nfilt / mapped_reads * 100`). |

## References
Reimplementation of the workflow used in https://doi.org/10.1016/j.watres.2023.119919.
