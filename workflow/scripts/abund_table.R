# see https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#r-and-r-markdown

#load required package
suppressPackageStartupMessages({
  if (!require("data.table")) {
    install.packages("data.table")
    require("data.table")
  }
  if (!require("dplyr")) {
    install.packages("dplyr")
    require("dplyr")
  }
})

main <- function(
  max_threads,
  mapping_overviews,
  totalreads_files,
  totalfilteredreads_files,
  database_tax,
  output,
  minalignlen,
  minIDfrac
) {
  #set max threads for data.table
  setDTthreads(as.integer(max_threads))

  #aggregate totalreads.csv files
  total_reads <- rbindlist(
    lapply(
      totalreads_files,
      fread,
      header = FALSE,
      sep = ",",
      col.names = c("barcode", "total_reads"),
      colClasses = c("character", "integer")
    )
  )

  #aggregate totalreads_filtered.csv files
  total_filtered_reads <- rbindlist(
    lapply(
      totalfilteredreads_files,
      fread,
      header = FALSE,
      sep = ",",
      col.names = c("barcode", "total_filtered_reads"),
      colClasses = c("character", "integer")
    )
  )

  #read taxonomy
  tax_db <- fread(
    database_tax,
    header = FALSE,
    sep = "\t",
    colClasses = "character",
    col.names = c("OTU", "tax")
  )

  if (length(mapping_overviews) == 0L) {
    stop("There are no mapping files to process. Exiting...", call. = FALSE)
  }

  mappings <- lapply(mapping_overviews, function(file) {
    mapping <- fread(
      file,
      header = FALSE,
      sep = " ",
      col.names = c(
        "readID",
        "SAMflag",
        "OTU",
        "querylen",
        "alnlen",
        "MapIDfrac",
        "NMtag",
        "alnscore",
        "MinimapID"
      )
    )
    mapping
  })
  names(mappings) <- gsub("/.*$", "", gsub("\\..*$", "",basename(mapping_overviews)))

  mappings <- rbindlist(
    mappings,
    idcol = "barcode"
  )

  # merge total reads per barcode
  mappings <- total_reads[mappings, on = "barcode"]
  mappings <- total_filtered_reads[mappings, on = "barcode"]

  # Calc ratio in pct between alignment length and query sequence length
  mappings[, alnfrac := round(alnlen / querylen, 3)]

  # calc number of "mappings" per barcode before filtering
  mappings[, mapped_reads := .N, by = barcode]

  # create a stats summary per barcode
  summary <- mappings[
    ,
    .(
      median_alnlen = median(alnlen),
      median_querylen = median(querylen),
      median_IDfrac = median(MapIDfrac),
      mapped_reads = mapped_reads[1],
      total_reads = total_reads[1],
      total_filtered_reads = total_filtered_reads[1]
    ),
    by = barcode
  ]
  summary[, pct_mapped := round(mapped_reads / total_filtered_reads * 100, 2)]

  # filter those with too short alignment or too low ID
  nmappingsbefore <- as.numeric(nrow(mappings))
  minalignlen <- as.integer(minalignlen)
  mappings <- mappings[alnlen >= minalignlen]
  if(!(minIDfrac <= 1 && minIDfrac >= 0))
    stop("Minimum mapping identity must be given as a fraction between 0.0 and 1.0")
  mappings <- mappings[MapIDfrac >= minIDfrac]
  nmappingsafter <- as.numeric(nrow(mappings))
  nfiltered <- nmappingsbefore - nmappingsafter

  # message
  message("Total (raw) reads across all barcodes: ", total_reads[, sum(total_reads)])
  message("Total (filtered) reads across all barcodes: ", total_filtered_reads[, sum(total_filtered_reads)])
  if (nfiltered == 0) {
    warning("0 mappings have been filtered by alignment length. This is highly suspicious unless there are very few reads. Manual inspection may be a good idea.")
  } else {
    message(
      paste0(
        nfiltered,
        " mappings (",
        round(nfiltered / total_filtered_reads[, sum(total_filtered_reads)] * 100, digits = 2L),
        "% of total) with too short alignment (<",
        minalignlen,
        "bp), and/or too many mismatches (<",
        minIDfrac * 100L,
        "% identity) have been filtered \nBefore: ",
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
  fwrite(
    summary,
    file.path(output, "summary.txt"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  #join taxonomy with mapping
  joined <- dplyr::left_join(mappings[, c("barcode", "OTU")], tax_db, by = "OTU")

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
    file.path(output, "otutable_mappedreads.tsv"),
    sep = "\t",
    col.names = TRUE,
    na = "NA",
    quote = FALSE
  )

  #extract abundances to be able to normalise to total reads
  abund <- otutable[, -c("OTU", ..tax_levels)]

  #make sure the order of barcodes is identical between otutable+reads
  total_filtered_reads <- total_filtered_reads[barcode %chin% colnames(abund)]
  total_filtered_reads[, barcode := factor(barcode, colnames(abund))]
  total_filtered_reads <- total_filtered_reads[order(barcode)]

  #normalise to total reads per barcode, in pct, rounded
  abund[] <- sweep(
    abund,
    2,
    total_filtered_reads[, total_filtered_reads],
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
    file.path(output, "otutable_normalised.tsv"),
    sep = "\t",
    col.names = TRUE,
    na = "NA",
    quote = FALSE
  )
}

main(
  max_threads = snakemake@threads,
  mapping_overviews = snakemake@input[["mapping_overviews"]],
  totalreads_files = snakemake@input[["totalreads_files"]],
  totalfilteredreads_files = snakemake@input[["totalfilteredreads_files"]],
  database_tax = snakemake@config[["db_tax"]],
  output = snakemake@config[["output_dir"]],
  minalignlen = snakemake@config[["minalignlen"]],
  minIDfrac = snakemake@config[["minIDfrac"]]
)

# for debugging:
# run_dir <- "./"
# max_threads <- 4
# mapping_overviews <- list.files(file.path(run_dir, "tmp/samples/"), pattern = "idmapped.txt", recursive = TRUE, full.names = TRUE)
# totalreads_files <- list.files(file.path(run_dir, "tmp/samples/"), pattern = "_totalreads.csv", recursive = TRUE, full.names = TRUE)
# database_tax <- "/databases/midas/MiDAS5.3_20240320/tax_complete_qiime.txt"
# output <- file.path(run_dir, "output")
# minalignlen <- 800
# minIDfrac <- 0.99
