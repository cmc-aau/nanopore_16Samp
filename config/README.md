# Configuration

 - `input_dir`: The input folder is expected to contain a subfolder for each sampleID/barcode, then all fastq files in each subfolder are concatenated and the folder name is used as sample ID downstream. This is usually the "fastq_pass" folder from nanopore sequencing and basecalling output 

 - `output_dir`: Output directory with the final results

 - `tmp_dir`: Directory for temporary files

 - `log_dir`: Directory for log files for all rules

 - `db_fasta`: Database for minimap2 in fasta file format. Each sequence header must only contain an ID matching the taxonomy file below, nothing else.

 - `db_tax`: 2-column TSV file with corresponding taxonomy for each sequence in the above fasta file. The first column is for the sequence IDs, the second a semi-colon separated taxonomy string.

 - `minalignlen`: Minimum alignment length for the mapping. Any alignments shorter than this threshold will be filtered

 - `minIDfrac`: Minimum identity threshold of each mapping (value must be between `0.0-1.0`). Any alignments with an identity lower than this threshold will be filtered.

 - `filtlong_args`: Arguments passed on to the `filtlong` command for pre-filtering reads
 
 - `max_threads`: Max number of threads for any rule

Have a look in the `.test` directory for example files.

## Database files
If you want to use the SILVA database, you can use [this R script](https://github.com/KasperSkytte/bioscripts/blob/main/extract_qiimetax.R) to create the two requires database files. For the MiDAS database the files can be downloaded directly from [download section](https://midasfieldguide.org/guide/downloads) on the website.
