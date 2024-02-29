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
  database_tax,
  output,
  minalignlen
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
      col.names = c("barcode", "reads"),
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
        "MapID",
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

  # Calc ratio in pct between query sequence length and alignment length
  mappings[, Qr := round(querylen / alnlen, 3)]

  # calc number of "mappings" per barcode before filtering
  mappings[, mapped_reads := .N, by = barcode]

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
  minalignlen <- as.integer(minalignlen)
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
  total_reads <- total_reads[barcode %chin% colnames(abund)]
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
  database_tax = snakemake@config[["db_tax"]],
  output = snakemake@config[["output_dir"]],
  minalignlen = snakemake@config[["minalignlen"]]
)
