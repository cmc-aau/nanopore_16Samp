# Full-length 16S nanopore workflow

Bash script to map reads against a taxonomic database and produce an abundance table compatible with ampvis2.

## Requirements

* minimap2
* samtools
* R (with data.table+dplyr installed)
* (filtlong)

You can also use the conda `environment.yml` file to create a conda environment.

## Usage

```
$ bash 16S_nanopore.sh -h
Script to process (demultiplexed) nanopore amplicon data for 16S community profiling.

version: 1.4.3
Options:
  -i    (required) Input fastq_pass folder with subfolders for each barcode.
  -o    (required) Output folder.
  -t    Max number of max_threads to use. (Default: all available except 2)
  -h    Display this help text and exit.
  -v    Print version and exit.

Additional options can be set by exporting environment variables before running the script:
  - database_fasta: Path to fasta file used for mapping
  - database_tax: Path to QIIME formatted taxonomy file matching the fasta file above
  - memlimit: Memory soft-limit for minimap2 in GB, fx "10g". Default: all free memory except 4GB
  - minalignlen: Minimum alignment length between query and reference sequence. Mappings with an alignment length equal to or shorter than this length will be not count in the abundance table. Default: 200
  - minfastqsize: Minimum fastq file size of all combined reads in bytes per barcode, otherwise skip processing the barcode, fx "100000"
  - minquality (currently not used): Remove reads with a quality of less than X % using filtlong
```

## Minimal example
```
$ export database_fasta=databases/MiDAS4.8.1_20210702/output/FLASVs.fa
$ export database_tax=databases/MiDAS4.8.1_20210702/output/tax_complete_qiime.txt
$ bash 16S_nanopore.sh -i 20220607_1300_MN34950_FAQ30081_bf6ff119/fastq_pass -o output_20220608 -t 10
#################################################
Host name: hurtigholger
Current user name: kapper
System time: 2023-02-20 13:33:42 (Europe/Copenhagen)
Script: /home/kapper/projects/FL16S_nanopore/16S_nanopore.sh
Script version: 1.4.3 (available at https://github.com/cmc-aau/FL16S_nanopore)
Current working directory: /home/kapper/projects/FL16S_nanopore
Input folder: /home/kapper/projects/FL16S_nanopore/data/20220315_1508_MN34021_FAQ30455_20b6f51c/fastq_pass
Output folder: /home/kapper/projects/FL16S_nanopore/output_20230220
Options set from environment variables:
  - database_fasta: databases/MiDAS4.8.1/FLASVs.fa
  - database_tax: databases/MiDAS4.8.1/tax_complete_qiime.txt
  - memlimit: 21g
  - minalignlen: 200
  - minfastqsize: 100000
  - minquality (currently not used): 80
Max. number of threads: 32
Log file: /home/kapper/projects/FL16S_nanopore/output_20230220/log_20230220_133342.txt
#################################################

 *** [2023-02-20 13:33:42] script message: Mapping all reads from each barcode to database
 *** [2023-02-20 13:33:42] script message: Processing barcode: barcode01
 *** [2023-02-20 13:33:42] script message:     (barcode: barcode01): Decompressing reads if gzip'ed and concatenating into a single file...
 *** [2023-02-20 13:33:45] script message:     (barcode: barcode01): Replacing sequence headers in fastq file with sha-256 checksums of themselves...
 *** [2023-02-20 13:33:46] script message:     (barcode: barcode01): Counting total number of reads...
 *** [2023-02-20 13:33:46] script message:     (barcode: barcode01): 44000 total reads
 *** [2023-02-20 13:33:46] script message:     (barcode: barcode01): Mapping against database and filtering output...
[M::mm_idx_gen::4.132*1.43] collected minimizers
[M::mm_idx_gen::4.418*2.13] sorted minimizers
[M::main::4.489*2.10] loaded/built the index for 90164 target sequence(s)
[M::mm_mapopt_update::4.526*2.09] mid_occ = 9765
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 90164
[M::mm_idx_stat::4.551*2.08] distinct minimizers: 953183 (40.87% are singletons); average occurrences: 23.271; average spacing: 5.511; total length: 122233683
[M::worker_pipeline::76.055*30.03] mapped 12879 sequences
[M::worker_pipeline::146.192*30.94] mapped 12892 sequences
[M::worker_pipeline::215.676*31.24] mapped 12866 sequences
[M::worker_pipeline::245.522*31.29] mapped 5363 sequences
[M::main] Version: 2.24-r1122
[M::main] CMD: minimap2 -ax map-ont --cap-sw-mem=21g -t 32 --secondary=no -K20M databases/MiDAS4.8.1/FLASVs.fa output_20230220/barcode01/barcode01_renamed.fastq
[M::main] Real time: 245.632 sec; CPU: 7682.347 sec; Peak RSS: 1.961 GB
 *** [2023-02-20 13:37:52] script message:     (barcode: barcode01): Creating mapping overview...
 *** [2023-02-20 13:38:14] script message: Generating abundance table
Total reads across all barcodes: 44000
118 mappings (across all barcodes) with too short alignment (min. alignment length: 200) have been filtered 
Before: 43977
After: 43859
Summary per barcode:
     barcode median_alnlen median_querylen mapped_reads total_reads pct_mapped nfilt pct_filt
1: barcode01          1472            1575        43977       44000      99.95   118     0.27
 *** [2023-02-20 13:38:15] script message: Done! Time elapsed: 00h:04m:33s

```
