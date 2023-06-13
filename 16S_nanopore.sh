#!/usr/bin/env bash
# This BASH script is based on the template from https://github.com/kasperskytte/bash_template
# License is MIT, which means the author takes no responsibilities, but you can use it for anything you want

##### settings start ######
#exit when a command fails (use "|| :" after a command to allow it to fail)
set -o errexit

# exit when a pipe fails
set -o pipefail

#disallow clobbering (overwriting) of files
#set -o noclobber

#print exactly what gets executed (useful for debugging)
#set -o xtrace

#set DK timezone if not set already
if [ -z "$(env | grep '^TZ=')" ]
then
  TZ="Europe/Copenhagen"
fi

version="1.4.4"

# Use all logical cores except 2 unless adjusted by user
max_threads=${max_threads:-$(($(nproc)-2))}

# Maximum memory in GB for minimap2. If unset default is max available minus 4GB
memlimit=${memlimit:-"$(($(free -g | awk '/^Mem:/{print $7}') - 4))g"}

# Path to fasta file used for mapping
database_fasta=${database_fasta:-"/databases/midas/MiDAS4.8.1_20210702/output/FLASVs.fa"}

# Path to QIIME formatted taxonomy file matching the fasta file above
database_tax=${database_tax:-"/databases/midas/MiDAS4.8.1_20210702/output/tax_complete_qiime.txt"}

# Remove reads with a quality of less than X %. If unset default is 80
minquality=${minquality:-80}

# Minimum fastq file size in bytes
minfastqsize=${minfastqsize:-100000}

# Minimum alignment length
minalignlen=${minalignlen:-200}

#log file name
logfilename=${logfilename:-"log_$(date '+%Y%m%d_%H%M%S').txt"}

##### settings end #####

#function to show a standard error message
usageError() {
  echo "Invalid usage: $1" 1>&2
  echo ""
  eval "bash $0 -h"
}

#function to check if a folder is present and empty
checkFolder() {
  #check user arguments
  if [ ! $# -eq 1 ]
  then
    echo "Error: function must be passed exactly 1 argument" >&2
    exit 1
  fi
  if [ -d "$1" ]
  then
      scriptMessage "A directory named '$1' already exists and is needed for this script to run. Please backup or delete the folder."
      scriptMessage "Exiting script."
      exit 1
  else
    mkdir -p "$1"
  fi
}

#function to add timestamps to progress messages
scriptMessage() {
  #check user arguments
  if [ ! $# -eq 1 ]
  then
    echo "Error: function must be passed exactly 1 argument" >&2
    exit 1
  fi
  echo " *** [$(date '+%Y-%m-%d %H:%M:%S')] script message: $1"
}

#function to check if executable(s) are available in $PATH
checkCommand() {
  args="$*"
  exit="no"
  for arg in $args
  do
    if [ -z "$(command -v "$arg")" ]
    then
      echo "${arg}: command not found"
      exit="yes"
    fi
  done
  if [ $exit == "yes" ]
  then
    exit 1
  fi
}

#check for all required commands before anything
checkCommand samtools R minimap2 perl cat awk sed #filtlong

#fetch and check options provided by user
while getopts ":hi:t:vo:" opt; do
case ${opt} in
  h )
    echo "Script to process (demultiplexed) nanopore amplicon data for 16S community profiling."
    echo
    echo "version: $version"
    echo "Options:"
    echo "  -i    (required) Input fastq_pass folder with subfolders for each barcode."
    echo "  -o    (required) Output folder."
    echo "  -t    Max number of max_threads to use. (Default: all available except 2)"
    echo "  -h    Display this help text and exit."
    echo "  -v    Print version and exit."
    echo
    echo "Additional options can be set by exporting environment variables before running the script:"
    echo "  - database_fasta: Path to fasta file used for mapping"
    echo "  - database_tax: Path to QIIME formatted taxonomy file matching the fasta file above"
    echo "  - memlimit: Memory soft-limit for minimap2 in GB, fx \"10g\". Default: all free memory except 4GB"
    echo "  - minalignlen: Minimum alignment length between query and reference sequence. Mappings with an alignment length equal to or shorter than this length will be not count in the abundance table. Default: 200"
    echo "  - minfastqsize: Minimum fastq file size of all combined reads in bytes per barcode, otherwise skip processing the barcode, fx \"100000\""
    echo "  - minquality (currently not used): Remove reads with a quality of less than X % using filtlong"
    exit 1
    ;;
  i )
    input=$OPTARG
    ;;
  t )
    max_threads=$OPTARG
    ;;
  v )
    echo "version: $version"
    exit 0
    ;;
  o )
    output=$OPTARG
    ;;
  \? )
    usageError "Invalid Option: -$OPTARG"
    exit 1
    ;;
  : )
    usageError "Option -$OPTARG requires an argument"
    exit 1
    ;;
esac
done
shift $((OPTIND -1)) #reset option pointer

#check required options
if [[ -z "$input" ]]
then
  usageError "option -i is required"
	exit 1
fi
if [[ -z "$output" ]]
then
	usageError "option -o is required"
	exit 1
fi
if ! [ "$minalignlen" -eq "$minalignlen" ] 2> /dev/null
then
  usageError "minalignlen must be an integer"
  exit 1
fi

# Start workflow

#check files first
if [ ! -s "$database_fasta" ]
then
  echo "Error: database fasta file doesn't exist or is empty"
  exit 1
fi

if [ ! -s "$database_tax" ]
then
  echo "Error: database taxonomy file doesn't exist or is empty"
  exit 1
fi

demuxfolders=$(
  find -L \
    "${input}"/ \
    -maxdepth 1 \
    -mindepth 1 \
    -type d \
    ! -iregex ".*unclassified$" \
    ! -iregex ".*unknown$"
)

if [ -z "$demuxfolders" ]
then
  echo "Error: There are no barcode subfolders in ${input}"
  exit 1
fi

checkFolder "${output}"

logfilepath="${output}/${logfilename}"


#wrap everything in a function to allow writing stderr+stdout to log file
main() {
  echo "#################################################"
  echo -n "Host name: $(hostname)"
  if [ -f "/.dockerenv" ]
  then
    echo -n " (Docker container)"
  fi
  echo
  echo "Current user name: $(whoami)"
  echo "System time: $(date '+%Y-%m-%d %H:%M:%S') (${TZ})"
  echo "Script: $(realpath "$0")"
  echo "Script version: ${version} (available at https://github.com/cmc-aau/FL16S_nanopore)"
  echo "Current working directory: $(pwd)"
  echo "Input folder: $(realpath "$input")"
  echo "Output folder: $(realpath "$output")"
  echo "Options set from environment variables:"
  echo "  - database_fasta: ${database_fasta}"
  echo "  - database_tax: ${database_tax}"
  echo "  - memlimit: ${memlimit}"
  echo "  - minalignlen: ${minalignlen}"
  echo "  - minfastqsize: ${minfastqsize}"
  echo "  - minquality (currently not used): ${minquality}"
  echo "Max. number of threads: ${max_threads}"
  echo "Log file: $(realpath "$logfilepath")"
  echo "#################################################"
  echo

  # File for summarizing number of reads per barcode
  total_reads_file="${output}/totalreads.csv"
  #clear total reads file
  true > "$total_reads_file"

  scriptMessage "Mapping all reads from each barcode to database"
  for barcode in $demuxfolders
  do
    barcodename=$(basename "$barcode")
    scriptMessage "Processing barcode: ${barcodename}"

    barcodefolder="${output}/${barcodename}"
    mkdir -p "$barcodefolder"

    fastqfiles=$(find -L "$barcode" -type f \( -iname '*.f*q' -o -iname '*.f*q.gz' \))

    #continue with next folder if no fastq files found
    if [ -z "$fastqfiles" ]
    then
      scriptMessage "    (barcode: ${barcodename}): No fastq files found, skipping..."
      continue
    fi

    #all fastq files will be concatenated into a single one
    barcode_allreads="${barcodefolder}/${barcodename}.fastq"
    barcode_allreads_renamed="${barcodefolder}/${barcodename}_renamed.fastq"
    #but first check if barcode_allreads_renamed.fastq already exists (empty or not)
    if [ -r "$barcode_allreads" ]
    then
      scriptMessage "    (barcode: ${barcodename}): barcode_allreads.fastq file already exists, skipping concatenation of all reads..."
      # #use find -mindepth 1 to avoid removing the barcode folder itself
      # #shellcheck disable=SC2046
      # rm -rf $(find "$barcodefolder" -mindepth 1 | grep -v 'barcode_allreads_renamed.fastq$') #don't quote
    elif [ -n "$fastqfiles" ]
    then
      #decompress if files are gzip'ed
      scriptMessage "    (barcode: ${barcodename}): Decompressing reads if gzip'ed and concatenating into a single file..."
      #dont quote $fastqfiles, it must be expanded
      #shellcheck disable=SC2086
      gunzip -cdfq $fastqfiles > "$barcode_allreads"

      scriptMessage "    (barcode: ${barcodename}): Replacing sequence headers in fastq file with sha-256 checksums of themselves..."
      seqid_rename_file="${barcodefolder}/seqid_checksums.tsv"
      #replaces fastq headers with a sha-256 checksum of them
      #to shorten them so samtools doesn't complain its too long
      #old headers and new headers are written to a CSV file to be able to backtrack
      awk 'NR%4 == 1 {print $0}' "$barcode_allreads" |\
        perl -MDigest::SHA=sha256_hex -nlE'say"$_\t".sha256_hex($_)' >\
        "$seqid_rename_file"

      awk 'FNR == NR {
        header[FNR]="@"$NF
        next
        }
        {
          if (FNR%4 == 1) {
            print header[(FNR-1)/4+1]
          }
          if (FNR%4 != 1) {
            print $0
          }
        }' \
        "$seqid_rename_file" \
        "$barcode_allreads" >\
        "$barcode_allreads_renamed"
    fi

    #first check if there is enough data
    filesize=$(stat -c%s "$barcode_allreads_renamed")
    if (( filesize < minfastqsize ))
    then
      scriptMessage "    (barcode: ${barcodename}): Insufficient or no data, skipping..."
      continue
    else
      mappingfile="${barcodefolder}/${barcodename}.idmapped.txt"

      scriptMessage "    (barcode: ${barcodename}): Counting total number of reads..."
      num_reads=$(($(wc -l < "$barcode_allreads_renamed") / 4))
      echo "${barcodename},${num_reads}" >> "$total_reads_file"
      scriptMessage "    (barcode: ${barcodename}): ${num_reads} total reads"

      if [ -s "$mappingfile" ]
      then
        scriptMessage "    (barcode: ${barcodename}): Mapping file has already been generated, skipping..."
      else
        # Map reads to reference database
        scriptMessage "    (barcode: ${barcodename}): Mapping against database and filtering output..."
        minimap2 \
          -ax map-ont \
          --cap-sw-mem="$memlimit" \
          -t "$max_threads" \
          --secondary=no \
          "$database_fasta" \
          -K20M "$barcode_allreads_renamed" \
          | samtools view \
            -F 256 \
            -F 4 \
            -F 2048 \
            --threads "${max_threads}" \
            -o "${barcodefolder}/${barcodename}.sam"

        # Create mapping overview
        scriptMessage "    (barcode: ${barcodename}): Creating mapping overview..."
        sed '/^@/ d' "${barcodefolder}/${barcodename}.sam" | \
          awk '{
            for(i=1;i<=NF;i++){
              if($i ~ /^NM:i:/){sub("NM:i:", "", $i); mm = $i}
            }
            split($6, count, /[^0-9]+/);
            split($6, type, /[^A-Z]*/);
            for(i=1; i<= length(count)-1; i++){
              if(type[i + 1] ~ /[DIM]/){aln+=count[i]};
            }
            print $1, $2, $3, length($10), aln, (aln - mm)/aln, $12, $14, $20
            aln=0;
          }' \
          > "$mappingfile"
      fi
    fi
  done

  # Filter reads based on quality
  #scriptMessage "Filtering low quality reads <${minquality}"
  #for f in $input/concatenated/*.fastq; do
  #barcodename=$(basename "$f" .fastq)
  #if [ -s output/filtered/$barcodename.filtered.fastq ]; then echo "output/filtered/$barcodename.filtered.fastq has already been generated";  
  #else
  #filtlong --min_mean_q $minquality $f > \
  #  output/filtered/$barcodename.filtered.fastq
  #fi

  scriptMessage "Generating abundance table"
  R --slave --args "${max_threads}" "${database_tax}" "${output}" "${total_reads_file}" "${minalignlen}" << 'makeOTUtable'
    #extract passed args from shell script
    args <- commandArgs(trailingOnly = TRUE)

    # for debugging if it fails
    # args <- list(
    #   max_threads = 32,
    #   database_tax = "databases/MiDAS4.8.1/tax_complete_qiime.txt",
    #   output = "/user_data/ksa/projects/midas_onlinekontrol/nanopore_runs/2023-may/p1/",
    #   total_reads_file = "/user_data/ksa/projects/midas_onlinekontrol/nanopore_runs/2023-may/p1/totalreads.csv",
    #   minalignlen = 200
    # )

    #load required package
    suppressPackageStartupMessages({
      if (!require("data.table")) {
        install.packages("data.table")
        require("data.table")
      }
    })

    #set max threads for data.table
    setDTthreads(as.integer(args[[1]]))

    #get total reads per sample
    total_reads <- fread(
      args[[4]],
      header = FALSE,
      sep = ",",
      col.names = c("barcode", "reads"),
      colClasses = c("character", "integer")
    )

    #read taxonomy
    tax_db <- fread(
      args[[2]],
      header = FALSE,
      sep = "\t",
      colClasses = "character",
      col.names = c("OTU", "tax")
    )

    files <- list.files(
      path = args[[3]],
      recursive = TRUE,
      pattern = ".idmapped.txt"
    )

    if (length(files) == 0L) {
      stop("There are no mapping files to process. Exiting...", call. = FALSE)
    }

    mappings <- lapply(files, function(file) {
      mapping <- fread(
        file.path(args[[3]], file),
        header = FALSE,
        sep = " ",
        col.names = c(
          "readID",
          "SAMflag",
          "OTU",
          "querylen",
          "alnlen",
          "MapID",
          "NMtag",
          "alnscore",
          "MinimapID"
        )
      )
      mapping
    })
    names(mappings) <- gsub("/.*$", "", files)

    mappings <- rbindlist(
      mappings,
      idcol = "barcode"
    )

    # merge total reads per barcode
    mappings <- mappings[total_reads, on = "barcode"]

    # Calc ratio in pct between query sequence length and alignment length
    mappings[, Qr := round(querylen / alnlen, 3)]

    # calc number of "mappings" per barcode before filtering
    mappings[, mapped_reads := .N, by = barcode]

    # Write out detailed mappings before filtering
    fwrite(
      mappings,
      file.path(args[[3]], "mappings_detailed.txt"),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )

    # create a stats summary per barcode
    summary <- mappings[
      ,
      .(
        median_alnlen = median(alnlen),
        median_querylen = median(querylen),
        mapped_reads = mapped_reads[1],
        total_reads = reads[1]
      ),
      by = barcode
    ]
    summary[, pct_mapped := round(mapped_reads / total_reads * 100, 2)]

    #filter those with too short alignment and message
    nmappingsbefore <- as.numeric(nrow(mappings))
    minalignlen <- as.integer(args[[5]])
    mappings <- mappings[alnlen >= minalignlen]
    nmappingsafter <- as.numeric(nrow(mappings))

    message("Total reads across all barcodes: ", total_reads[, sum(reads)])

    if (nmappingsbefore == nmappingsafter) {
      message(
        "0 mappings have been filtered (min. alignment length: ",
        minalignlen,
        ")."
      )
    } else {
      message(
        paste0(
          nmappingsbefore - nmappingsafter,
          " mappings (across all barcodes) with too short alignment (min. alignment length: ",
          minalignlen,
          ") have been filtered \nBefore: ",
          nmappingsbefore,
          "\nAfter: ",
          nmappingsafter
        )
      )
    }

    # calc number of "mappings" per barcode after filtering and print+write out summary
    mappings[, nfilt := mapped_reads - .N, by = barcode]
    summary <- summary[unique(mappings[, c("barcode", "nfilt")]), on = "barcode"]
    summary[, pct_filt := round(nfilt / mapped_reads * 100, 2)]
    message("Summary per barcode:")
    print(summary)
    fwrite(
      summary,
      file.path(args[[3]], "summary.txt"),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )

    #join taxonomy with mapping
    joined <- tax_db[mappings[, c("barcode", "OTU")], on = "OTU"]

    #transform into a merged "OTU table",
    #includes both abundances and taxonomy (old school ampvis format)
    otutable <- dcast(
      joined,
      OTU + tax ~ barcode,
      fun.aggregate = length
    )

    tax_levels <- c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )

    #split up taxonomy column into separate columns with tax levels
    otutable[, c(tax_levels) := tstrsplit(tax, ";", fixed = FALSE)]
    otutable[, tax := NULL]

    #write out non-normalised table
    fwrite(
      otutable,
      file.path(args[[3]], "otutable_mappedreads.tsv"),
      sep = "\t",
      col.names = TRUE,
      na = "NA",
      quote = FALSE
    )

    #extract abundances to be able to normalise to total reads
    abund <- otutable[, -c("OTU", ..tax_levels)]

    #make sure the order of barcodes is identical between otutable+reads
    total_reads[, barcode := factor(barcode, colnames(abund))]
    total_reads <- total_reads[order(barcode)]

    #normalise to total reads per barcode, in pct, rounded
    abund[] <- sweep(
      abund,
      2,
      total_reads[, reads],
      "/"
    ) * 100
    abund[] <- round(abund, digits = 4)

    #stitch together normalised otutable and write out
    otutable_norm <- data.table(
      otutable[, "OTU"],
      abund,
      otutable[, ..tax_levels]
    )
    fwrite(
      otutable_norm,
      file.path(args[[3]], "otutable_normalised.tsv"),
      sep = "\t",
      col.names = TRUE,
      na = "NA",
      quote = FALSE
    )
makeOTUtable

  duration=$(printf '%02dh:%02dm:%02ds\n' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60)))
  scriptMessage "Done! Time elapsed: $duration"
}

#clear log file first if it happens to already exist
true > "$logfilepath"
main |& tee -a "$logfilepath"
