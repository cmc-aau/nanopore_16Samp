#!/usr/bin/env bash

##### settings start ######
#exit when a command fails (use "|| :" after a command to allow it to fail)
set -o errexit

# exit when a pipe fails
set -o pipefail

#disallow undeclared variables
set -o nounset

#disallow clobbering (overwriting) of files
#set -o noclobber

#print exactly what gets executed (useful for debugging)
#set -o xtrace

version="1.0"

#use all logical cores except 2 unless adjusted by user
max_threads=${max_threads:-$(($(nproc)-2))}

#maximum memory for minimap2. If unset default is max available minus 4GB
memlimit=${memlimit:-$(($(free -g | awk '/^Mem:/{print $7}') - 4))}

# Path to fasta file used for mapping
database_fasta=${database_fasta:-"/space/databases/midas/MiDAS4.8.1_20210702/output/FLASVs.fa"}

# Path to QIIME formatted taxonomy file matching the fasta file above
database_tax=${database_tax:-"/space/databases/midas/MiDAS4.8.1_20210702/output/tax_complete_qiime.txt"}

# Remove reads with a quality of less than X %. If unset default is 80
minquality=${minquality:-80}

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
  if [ -d $1 ]
  then
      scriptMessage "A directory named '$1' already exists and is needed for this script to run. Please backup or delete the folder."
      scriptMessage "Exiting script."
      exit 1
  else
    mkdir -p $1
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
  args="$@"
  exit="no"
  for arg in $args
  do
    if [ -z $(command -v $arg) ]
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
#flags for required options, checked after getopts loop
i_flag=0
o_flag=0
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
    i_flag=1
    ;;
  t )
    MAX_max_threads=$OPTARG
    ;;
  v )
    echo "version: $version"
    exit 0
    ;;
  o )
    output=$OPTARG
    o_flag=1
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

#check all required options
if [ $i_flag -eq 0 ]
then
	usageError "option -i is required"
	exit 1
fi
if [ $o_flag -eq 0 ]
then
	usageError "option -o is required"
	exit 1
fi

##### START OF ACTUAL SCRIPT #####
#print settings here?
#scriptMessage "Running Nanopore raw mapping script for demultiplexed files (max max_threads $max_threads)"

# Start workflow
database_name=$(echo "$database_fasta" | grep -o "[^/]*$" | grep -o "^[^\.]*")

#decompress if files are gzip'ed
find $input -iname '*.gz' -exec gunzip -q {} \;

fastqfiles=$(find $input -iname '*.f*q' -printf '%h\n' | sort -u)
if [ -z "$fastqfiles" ]
then
  echo "Error: no fastq files found in ${input}, exiting..."
  exit 1
fi

checkFolder output/
checkFolder output/mapped
checkFolder output/filtered/

for dirpath in $fastqfiles/
do
  barcodename=$(basename "$dirpath")
  if [ -s $output/concatenated/$barcodename.fastq ]
  then
    echo "$barcodename has already been concatenated";  
  else
    (
    mkdir -p $output/concatenated/
    echo "Concatenating $barcodename"
      cat $dirpath/*.f*q | \
      sed 's/\(runid.*start_\)//' | \
      tr -d '[:blank:]' > $output/concatenated/$barcodename.fastq
    )
  fi
done  

# Filter reads based on quality
#for f in $input/concatenated/*.fastq; do
#filename=$(basename "$f" .fastq)
#if [ -s output/filtered/$filename.filtered.fastq ]; then echo "output/filtered/$filename.filtered.fastq has already been generated";  
#else
#filtlong --min_mean_q $minquality $f > \
#  output/filtered/$filename.filtered.fastq
#fi


for file in $output/concatenated/*.fastq; do
filename=$(basename "$file" .fastq)
if [ -s output/$filename.idmapped.txt ]
then
  echo "output/$filename.idmapped.txt has already been generated"
else
  # Map reads to reference database
  minimap2 -ax map-ont \
    --cap-sw-mem=$memlimit \
    -t $max_threads \
    --secondary=no \
    $database_fasta \
    -K20M $file > \
    output/mapped/$filename.sam

  samtools view -F 256 \
    -F 4 \
    -F 2048 output/mapped/$filename.sam \
    -o output/mapped/$filename"_nodupes.sam"
  # Create mapping overview
  scriptMessage "Creating mapping overview"
  sed '/^@/ d' output/mapped/${filename}_nodupes.sam | \
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
    sed 's/time=/\t/' | sed 's/Zflow_cell.*barcode=/\t/' > output/$filename.idmapped.txt
fi
done

scriptMessage "Generating OTU/mapping table for loading into ampvis2"

R --slave --args "$max_threads" "$database_tax" "$database_name" << 'makeOTUtable'
  #extract passed args from shell script
  args <- commandArgs(trailingOnly = TRUE)
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

  setDTmax_threads(as.integer(args[[1]]))

  #read taxonomy
  taxDB <- fread(args[[2]],
                header = FALSE,
                sep = "\t",
                colClasses = "character",
                col.names = c("OTU", "tax"))

  #read mappings
  mappings<-data.frame(OTU=NULL,SeqID=NULL)
  files<-list.files(path = "output/", pattern = ".idmapped.txt")
  for (file in files){
    mapping <- fread(
      paste0("output/",file),
      header = FALSE,
      #select = c(3),
      #colClasses = "character",
    col.names = c("readID", "read_time", "add_info")
    ) %>% separate(add_info, c("barcode", "SAMflag", "OTU", "Qlen", "alnlen", "MapID", "NMtag", "alnscore", "MinimapID"), " ") %>% mutate(SeqID=sub(".idmapped.txt","",file))
    mappings<-rbind(mappings,mapping)
  }

  # Subset mappings based on 
  mappings_s <- mappings %>% mutate(Qr = as.numeric(Qlen)/as.numeric(alnlen)) %>% 
  subset(., Qr < 1.15 & Qr > 0.85)
  mappings_to_otu <- mappings_s %>% 
  select(c("SeqID", "OTU"))

  # Write mappings out for further analysis in R
  question <- askYesNo("Do you want to write out a detailed mapping file?", default=TRUE, prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
  if (question == TRUE) {
    fwrite(mappings_s, "output/mappings_filt15.txt", quote = F, sep = "\t", row.names = F, col.names = T)
  }else if (question == FALSE) {
    print("Okay, we'll resume..")
  }

  # Define function for 
  #join taxonomy with mapping
  joined <- taxDB[mappings_to_otu, on = "OTU"]

  #transform into "OTU table", where rows are OTU's, columns are sample AND taxonomy
  BIOMotutable <- dcast(joined, OTU + tax ~ SeqID, fun.aggregate = length) %>% setDT()
  BIOMotutable[,taxonomy := gsub("[a-z]__", "", tax)]
  BIOMotutable[,tax := NULL]
  BIOMotutable[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") := tstrsplit(taxonomy, ";", fixed=FALSE)]
  BIOMotutable[,taxonomy := NULL]
  #write out
  fwrite(BIOMotutable, paste0("output/otutable_", args[[3]], ".txt"), sep = "\t", col.names = TRUE, na = "NA", quote = FALSE)
makeOTUtable

duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
scriptMessage "Done in: $duration, enjoy!"

##### END OF ACTUAL SCRIPT #####

#print elapsed time since script was invoked
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
scriptMessage "Done in: $duration!"
exit 0