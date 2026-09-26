#!/usr/bin/env Rscript

# ===============================================================================
# spikein_barcode_merge.R
#
# Created for the grothlab/crepas pipeline by:
#     - Samuel Ruiz-Perez <samper@cancer.dk>
#     https://github.com/grothlab/crepas/
#
# Description:
#     Merge one or more spike-in barcode count tables by summing R1_count and
#     R2_count for each unique barcode_target/barcode_id/barcode_sequence.
# ===============================================================================

options(show.error.locations = TRUE)

required.libs <- c("readr", "dplyr", "argparse")

for (lib in required.libs) {
  suppressPackageStartupMessages({
    if (!require(lib, character.only = TRUE)) {
      stop("Missing required R package: ", lib)
    }
  })
}

# ===============================================================================
# Argument parsing
# ===============================================================================

parser <- argparse::ArgumentParser(
  description = "Merge spike-in barcode count TSVs across technical replicates"
)

parser$add_argument(
  "-c", "--count_tables",
  action = "store",
  type = "character",
  nargs = "+",
  required = TRUE,
  help = "One or more count tables with columns: barcode_target, barcode_id, barcode_sequence, R1_count, R2_count"
)

parser$add_argument(
  "-p", "--prefix",
  action = "store",
  type = "character",
  default = "spikein_barcodes",
  help = "Prefix for output file [default: spikein_barcodes]"
)

parser$add_argument(
  "-o", "--outdir",
  action = "store",
  type = "character",
  default = "./",
  help = "Output directory for merged count table [default: ./]"
)

opt <- parser$parse_args()
opt_count_tables <- opt$count_tables
opt_prefix <- opt$prefix
opt_outdir <- opt$outdir

# ===============================================================================
# Helpers
# ===============================================================================

required_cols <- c(
  "barcode_target",
  "barcode_id",
  "barcode_sequence",
  "R1_count",
  "R2_count"
)

read_count_table <- function(path) {
  if (!file.exists(path)) {
    stop("Count table not found: ", path)
  }

  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Table '", path, "' is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df
}

# ===============================================================================
# Merge tables
# ===============================================================================

message("[", Sys.time(), "] Reading ", length(opt_count_tables), " count table(s)...")

table_list <- lapply(opt_count_tables, read_count_table)
all_counts <- dplyr::bind_rows(table_list)

merged_counts <- all_counts %>%
  dplyr::group_by(barcode_target, barcode_id, barcode_sequence) %>%
  dplyr::summarise(
    R1_count = sum(R1_count, na.rm = TRUE),
    R2_count = sum(R2_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(barcode_target, barcode_id)

if (!dir.exists(opt_outdir)) {
  dir.create(opt_outdir, recursive = TRUE)
}

output_file <- file.path(opt_outdir, paste0(opt_prefix, ".tsv"))
readr::write_tsv(merged_counts, output_file, col_names = TRUE)

message("[", Sys.time(), "] Wrote merged count table: ", output_file)
message("[", Sys.time(), "] Rows: ", nrow(merged_counts))
