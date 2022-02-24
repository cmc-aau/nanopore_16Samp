#!/usr/bin/env bash

##### settings start ######
#exit when a command fails (use "|| :" after a command to allow it to fail)
set -o errexit

# exit when a pipe fails
set -o pipefail

#disallow clobbering (overwriting) of files
#set -o noclobber

#print exactly what gets executed (useful for debugging)
#set -o xtrace

version="1.2"

#use all logical cores except 2 unless adjusted by user
max_threads=${max_threads:-$(($(nproc)-2))}

#maximum memory in GB for minimap2. If unset default is max available minus 4GB
memlimit=${memlimit:-"$(($(free -g | awk '/^Mem:/{print $7}') - 4))g"}

# Path to fasta file used for mapping
database_fasta=${database_fasta:-"/space/databases/midas/MiDAS4.8.1_20210702/output/FLASVs.fa"}

# Path to QIIME formatted taxonomy file matching the fasta file above
database_tax=${database_tax:-"/space/databases/midas/MiDAS4.8.1_20210702/output/tax_complete_qiime.txt"}

# Remove reads with a quality of less than X %. If unset default is 80
minquality=${minquality:-80}

# minimum fastq file size in bytes
minfastqsize=${minfastqsize:-100000}

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
checkCommand samtools R minimap2 filtlong

#fetch and check options provided by user
while getopts ":hi:t:vo:" opt; do
case ${opt} in
  h )
    echo "Some help text"
    echo "version: $version"
    echo "Options:"
    echo "  -h    Display this help text and exit."
    echo "  -v    Print version and exit."
    echo "  -i    (required) Input fastq_pass folder."
    echo "  -o    (required) Output folder."
    echo "  -t    Max number of max_threads to use. (Default: all available except 2)"
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

# Start workflow

# generate a name for the database based on its file name
database_name=$(
  echo "${database_fasta}" |\
    grep -o "[^/]*$" |\
    grep -o "^[^\.]*"
)

checkFolder "${output}"

#decompress if files are gzip'ed
scriptMessage "Decompressing f*q files if gzip'ed"
find "${input}" -iname '*.f*q.gz' -exec gunzip -q {} \;

demuxfolders=$(
  find "${input}"/* \
    -maxdepth 0 \
    -type d \
    ! -iregex ".*unclassified$" \
    ! -iregex ".*unknown$"
)

scriptMessage "Mapping all reads from each barcode to database"
for barcode in $demuxfolders
do
  barcodename=$(basename "$barcode")
	scriptMessage "Processing barcode: ${barcodename}"

  barcodefolder="${output}/${barcodename}"
  mkdir -p "$barcodefolder"

  fastqfiles=$(find "$barcode" -type f -iname '*.f*q')

	#continue with next folder if no fastq files found
	if [ -z "$fastqfiles" ]
	then
	  scriptMessage "    (barcode: ${barcodename}): No fastq files found, skipping..."
	  continue
	fi

	#concatenate all fastq files into one
	barcode_allreads="${barcodefolder}/${barcodename}.fastq"
	#but first check if barcode_allreads.fastq already exists (empty or not)
	if [ -r "$barcode_allreads" ]
	then
	  scriptMessage "    (barcode: ${barcodename}): barcode_allreads.fastq file already exists, skipping concatenation of all reads..."
	  # #use find -mindepth 1 to avoid removing the barcode folder itself
		# #shellcheck disable=SC2046
	  # rm -rf $(find "$barcodefolder" -mindepth 1 | grep -v 'barcode_allreads.fastq$') #don't quote
  elif [ -n "$fastqfiles" ]
	then
    scriptMessage "    (barcode: ${barcodename}): Concatenating all fastq files into a single file..."
		#dont quote $fastqfiles
    #shellcheck disable=SC2086
    cat $fastqfiles |\
      sed 's/\(runid.*start_\)//' |\
      tr -d '[:blank:]' > "$barcode_allreads"
	fi

	#first check if there is enough data
	filesize=$(stat -c%s "$barcode_allreads")
	if (( filesize < minfastqsize ))
	then
    scriptMessage "    (barcode: ${barcodename}): Insufficient or no data, skipping..."
    continue
  else
    mappingfile="${barcodefolder}/${barcodename}.idmapped.txt"
    if [ -s "$mappingfile" ]
    then
      scriptMessage "    (barcode: ${barcodename}): Mapping file has already been generated, skipping..."
    else
      # Map reads to reference database
      scriptMessage "    (barcode: ${barcodename}): Mapping against database..."
      minimap2 \
        -ax map-ont \
        --cap-sw-mem="$memlimit" \
        -t "$max_threads" \
        --secondary=no \
        "$database_fasta" \
        -K20M "$barcode_allreads" > \
        "${barcodefolder}/${barcodename}.sam"
      
      scriptMessage "    (barcode: ${barcodename}): Filtering mapping output..."
      samtools view \
        -F 256 \
        -F 4 \
        -F 2048 "${barcodefolder}/${barcodename}.sam" \
        --threads "${max_threads}" \
        -o "${barcodefolder}/${barcodename}_nodupes.sam"

      # Create mapping overview
      scriptMessage "    (barcode: ${barcodename}): Creating mapping overview..."
      sed '/^@/ d' "${barcodefolder}/${barcodename}_nodupes.sam" | \
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
        }' | \
        sed 's/time=/\t/' |\
        sed 's/Zflow_cell.*barcode=/\t/' \
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
R --slave --args "${max_threads}" "${database_tax}" "${database_name}" "${output}" << 'makeOTUtable'
  #extract passed args from shell script
  args <- commandArgs(trailingOnly = TRUE)

  #load required package
  suppressPackageStartupMessages({
    if(!require("data.table")) {
      install.packages("data.table")
      require("data.table")
    }
    if(!require("dplyr")) {
      install.packages("dplyr")
      require("dplyr")
    }
    if(!require("tidyr")) {
      install.packages("tidyr")
      require("tidyr")
    }
  })

  #set max threads for data.table
  setDTthreads(as.integer(args[[1]]))

  #read taxonomy
  taxDB <- fread(
    args[[2]],
    header = FALSE,
    sep = "\t",
    colClasses = "character",
    col.names = c("OTU", "tax")
  )

  #read mappings
  mappings <- data.frame(OTU=NULL, SeqID=NULL)

  files <- list.files(
    path = args[[4]],
    recursive = TRUE,
    pattern = ".idmapped.txt"
  )

  for (file in files) {
    mapping <- fread(
      paste0(args[[4]], "/", file),
      header = FALSE,
      #select = c(3),
      #colClasses = "character",
      col.names = c("readID", "read_time", "add_info")
    ) %>% 
      separate(
        add_info,
        c(
          "barcode",
          "SAMflag",
          "OTU",
          "Qlen",
          "alnlen",
          "MapID",
          "NMtag",
          "alnscore",
          "MinimapID"
        ),
        " ") %>%
      mutate(SeqID = sub(".idmapped.txt", "", file))
    mappings <- rbind(mappings, mapping)
  }

  # Subset mappings based on 
  mappings_s <- mappings %>%
    mutate(Qr = as.numeric(Qlen) / as.numeric(alnlen)) %>%
    subset(Qr < 1.15 & Qr > 0.85)
  mappings_to_otu <- mappings_s %>%
    select(c("SeqID", "OTU"))

  # Write out detailed mappings
  fwrite(
    mappings_s,
    paste0(args[[4]], "/mappings_detailed.txt"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Define function for 
  #join taxonomy with mapping
  joined <- taxDB[mappings_to_otu, on = "OTU"]

  #transform into "OTU table", where rows are OTU's, columns are sample AND taxonomy
  BIOMotutable <- dcast(
    joined,
    OTU + tax ~ SeqID,
    fun.aggregate = length
  )
  setDT(BIOMotutable)
  BIOMotutable[,taxonomy := gsub("[a-z]__", "", tax)]
  BIOMotutable[,tax := NULL]
  BIOMotutable[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") := tstrsplit(taxonomy, ";", fixed=FALSE)]
  BIOMotutable[,taxonomy := NULL]
  #write out
  fwrite(
    BIOMotutable,
    paste0(args[[4]], "/otutable_", args[[3]], ".txt"),
    sep = "\t",
    col.names = TRUE,
    na = "NA",
    quote = FALSE
  )
makeOTUtable

duration=$(printf '%02dh:%02dm:%02ds\n' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60)))
scriptMessage "Done! Time elapsed: $duration"
