# Full-length 16S nanopore workflow
Bash script to map reads against a taxonomic database and produce a abundance table compatible with ampvis2.

## Requirements
 - minimap2
 - samtools
 - R
 - (filtlong)

## Usage
```
$ bash workflow.sh -h
Some help text
version: 1.0
Options:
  -h    Display this help text and exit.
  -v    Print version and exit.
  -i    (required) Input fastq_pass folder.
  -o    (required) Output folder.
  -t    Max number of max_threads to use. (Default: all available except 2)
```
