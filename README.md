# Full-length 16S nanopore workflow

Bash script to map reads against a taxonomic database and produce an abundance table compatible with ampvis2.

## Requirements

* minimap2
* samtools
* R
* (filtlong)

On AAU just load the following modules before running the script:

```
module load Minimap2/2.15-foss-2018a
module load SAMtools/1.10-foss-2018a
module load R/3.5.0-foss-2018a-X11-20180131
```

## Usage

```
$ bash 16S_nanopore.sh -h
Script to process (demultiplexed) nanopore amplicon data for 16S community profiling.

version: 1.3.6
Options:
  -i    (required) Input fastq_pass folder with subfolders for each barcode.
  -o    (required) Output folder.
  -t    Max number of max_threads to use. (Default: all available except 2)
  -h    Display this help text and exit.
  -v    Print version and exit.

Additional options can be set by exporting environment variables before running the script:
  - database_fasta: Path to fasta file used for mapping
  - database_tax: Path to QIIME formatted taxonomy file matching the fasta file above
  - memlimit: Memory soft-limit for minimap2 in GB, fx "10g"
  - minfastqsize: Minimum fastq file size of all combined reads in bytes per barcode, otherwise skip processing the barcode, fx "100000"
  - minquality (currently not used): Remove reads with a quality of less than X % using filtlong

```

## Minimal example
```
$ export database_fasta=databases/MiDAS4.8.1_20210702/output/FLASVs.fa
$ export database_tax=databases/MiDAS4.8.1_20210702/output/tax_complete_qiime.txt
$ bash 16S_nanopore.sh -i 20220607_1300_MN34950_FAQ30081_bf6ff119/fastq_pass -o output_20220608 -t 10
 *** [2022-06-08 09:44:43] script message: Mapping all reads from each barcode to database
 *** [2022-06-08 09:44:43] script message: Processing barcode: barcode01
 *** [2022-06-08 09:44:43] script message:     (barcode: barcode01): Decompressing reads if gzip'ed and concatenating into a single file...
 *** [2022-06-08 09:44:43] script message:     (barcode: barcode01): Replacing sequence headers in fastq file with sha-256 checksums of themselves...
 *** [2022-06-08 09:44:43] script message:     (barcode: barcode01): Counting total number of reads...
 *** [2022-06-08 09:44:43] script message:     (barcode: barcode01): 5876 total reads
 *** [2022-06-08 09:44:43] script message:     (barcode: barcode01): Mapping against database...
[M::mm_idx_gen::2.328*1.32] collected minimizers
[M::mm_idx_gen::2.824*2.03] sorted minimizers
[M::main::2.834*2.03] loaded/built the index for 90164 target sequence(s)
[M::mm_mapopt_update::2.862*2.02] mid_occ = 9765
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 90164
[M::mm_idx_stat::2.884*2.01] distinct minimizers: 953183 (40.87% are singletons); average occurrences: 23.271; average spacing: 5.511; total length: 122233683
[M::worker_pipeline::59.005*9.34] mapped 5876 sequences
[M::main] Version: 2.24-r1122
[M::main] CMD: minimap2 -ax map-ont --cap-sw-mem=10g -t 10 --secondary=no -K20M databases/MiDAS4.8.1_20210702/output/FLASVs.fa output_20220608/barcode01/barcode01_renamed.fastq
[M::main] Real time: 59.056 sec; CPU: 551.376 sec; Peak RSS: 1.007 GB
 *** [2022-06-08 09:45:42] script message:     (barcode: barcode01): Filtering mapping output...
 *** [2022-06-08 09:45:43] script message:     (barcode: barcode01): Creating mapping overview...
 *** [2022-06-08 09:45:43] script message: Generating abundance table
 *** [2022-06-08 09:45:43] script message: Done! Time elapsed: 00h:01m:00s
```