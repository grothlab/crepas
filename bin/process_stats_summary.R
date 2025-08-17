#!/usr/bin/env Rscript
# Author: Samuel Ruiz-Pérez
# Script to process samtools_stats_summary table outputted by the grothlab/glseq pipeline

options(show.error.locations = TRUE)

required.libs <- c("dplyr", "readr", "tidyr", "ggplot2", "forcats", "argparse")

suppressPackageStartupMessages({
  lapply(required.libs, FUN = function(x) {
    do.call("require", list(x))
  })
})


parser <- ArgumentParser()


parser$add_argument("-s","--summary_table", action = "store",
                    type = "character",
                    help = "Summary table from samtools_stats_summary [required]")

parser$add_argument("-g","--endogenous_genome_name", action = "store",
                    type = "character",
                    default = "GRCh38",
                    help = "Name of the endogenous genome [required]")

parser$add_argument("-e","--exogenous_genome_name", action = "store",
                    type = "character",
                    default = "dm6",
                    help = "Name of the exogenous (spike-in) genome if applicable [optional]")

parser$add_argument("-n", "--prefix", action = "store",
                    default = "final_samtools_stats_summary",
                    type = "character",
                    help = "Prefix for output files and plot title")

parser$add_argument("-o", "--outdir", action = "store",
                    default = "./",
                    type = "character",
                    help = "Path to output directory for processed tables and plots")

opt <- parser$parse_args()

if (is.null(opt$summary_table)) {
  stop("The summary table is required.")
} else if (!file.exists(opt$summary_table)) {
  stop("The summary table file does not exist.")
}

if (is.null(opt$endogenous_genome_name)) {
  stop("The endogenous genome name is required.")
}

summary_table <- read_tsv(opt$summary_table, col_names = TRUE)

if (is.null(opt$exogenous_genome_name)) {
    
  tmp <- summary_table %>%
                  # In the column names replace any space or dash with underscore
                  rename_with(~ gsub("[- ]", "_", .x)) %>%
                  # Remove any row before library merging (containing ".Lb.")
                  filter(!grepl("\\.Lb\\.", ID)) %>%
                  # remove ".sorted" and ".sorted.bam" from the ID column
                  mutate(
                    ID = gsub("\\.sorted\\.bam$", "", ID),
                    ID = gsub("\\.sorted$", "", ID),
                    # The following is to avoid endo and exogenous files merging together
                    # ID = gsub("\\.flT2$", "_flT2", ID),
                    # ID = gsub("\\.allo$", "_allo", ID),
                    # ID = gsub("\\.flT3$", "_flT3", ID),
                    # ID = gsub("\\..shifted$", "_shifted", ID),
                    # ID = gsub("\\.dSp$", "_dSp", ID),
                    # split ID by ".", processing_step is the last element
                    processing_step = gsub(".*\\.", "", ID),
                    processing_step = fct_relevel(processing_step, rev(unique(processing_step))),
                    sample = gsub("\\..*$", "", ID)) %>%
                  select(-ID) %>%
                  # First, pivot the data to long format for easier reshaping
                  pivot_longer(
                    cols = -c(sample, processing_step),
                    names_to = "stat",
                    values_to = "value"
                  ) %>%
                  # sort by processing_step factor levels
                  arrange(processing_step) %>%
                  # Create new column names by combining processing_step and stat
                  mutate(stat = paste(processing_step, stat, sep = "_")) %>%
                  select(-processing_step) %>%
                  # Pivot back to wide format
                  pivot_wider(
                    names_from = stat,
                    values_from = value
                  )

  # If there is flT3_raw_total_sequences column, use it to calculate flT1_multimapping_reads
  if (any(grepl("flT3_raw_total_sequences", colnames(tmp)))) {
    tmp <- tmp %>%
          mutate(
            # Add a new column for flT1_multimapping_reads
            flT1_multimapping_reads = flT1_raw_total_sequences - flT3_raw_total_sequences) %>%
          # move flT1_multimapping_reads after last column with flT1 in the name
          relocate(flT1_multimapping_reads, .after = tail(grep("flT1", colnames(.[-which(colnames(.) == "flT1_multimapping_reads")])), n = 1))
  }                

} else {

  tmp <- summary_table %>%
                  # In the column names replace any space or dash with underscore
                  rename_with(~ gsub("[- ]", "_", .x)) %>%
                  # Remove any row before library merging (containing ".Lb.")
                  filter(!grepl("\\.Lb\\.", ID)) %>%
                  # remove ".sorted" and ".sorted.bam" from the ID column
                  mutate(
                    ID = gsub("\\.sorted\\.bam$", "", ID),
                    ID = gsub("\\.sorted$", "", ID),
                    # The following is to avoid endo and exogenous files merging together
                    ID = ifelse(grepl(paste0(".", opt$endogenous_genome_name, "."), ID),
                                gsub("\\.flT2$", paste0(".", opt$endogenous_genome_name, "_flT2"), ID), ID),
                    ID = ifelse(grepl(paste0(".", opt$exogenous_genome_name, "."), ID),
                                gsub("\\.flT2$", paste0(".", opt$exogenous_genome_name, "_flT2"), ID), ID),
                    ID = ifelse(grepl(paste0(".", opt$endogenous_genome_name, "."), ID),
                                gsub("\\.allo$", paste0(".", opt$endogenous_genome_name, "_allo"), ID), ID),
                    ID = ifelse(grepl(paste0(".", opt$exogenous_genome_name, "."), ID),
                                gsub("\\.allo$", paste0(".", opt$exogenous_genome_name, "_allo"), ID), ID),
                    ID = ifelse(grepl(paste0(".", opt$endogenous_genome_name, "."), ID),
                                gsub("\\..shifted$", paste0(".", opt$endogenous_genome_name, "_shifted"), ID), ID),
                    ID = ifelse(grepl(paste0(".", opt$exogenous_genome_name, "."), ID), 
                                gsub("\\..shifted$", paste0(".", opt$exogenous_genome_name, "_shifted"), ID), ID),
                    ID = ifelse(grepl(paste0(".", opt$endogenous_genome_name, "."), ID),
                                gsub("\\.flT3$", paste0(".", opt$endogenous_genome_name, "_flT3"), ID), ID),
                    ID = ifelse(grepl(paste0(".", opt$exogenous_genome_name, "."), ID),
                                gsub("\\.flT3$", paste0(".", opt$exogenous_genome_name, "_flT3"), ID), ID),
                    ID = ifelse(grepl(paste0("\\.", opt$endogenous_genome_name, "\\."), ID),
                          gsub("\\.dSp(.*)$", paste0(".", opt$endogenous_genome_name, "_dSp\\1"), ID), ID),
                    ID = ifelse(grepl(paste0("\\.", opt$exogenous_genome_name, "\\."), ID),
                          gsub("\\.dSp(.*)$", paste0(".", opt$exogenous_genome_name, "_dSp\\1"), ID), ID),
                    # split ID by ".", processing_step is the last element
                    processing_step = gsub(".*\\.", "", ID),
                    processing_step = fct_relevel(processing_step, rev(unique(processing_step))),
                    sample = gsub("\\..*$", "", ID)) %>%
                    dplyr::select(-ID) %>%
                    # First, pivot the data to long format for easier reshaping
                    pivot_longer(
                      cols = -c(sample, processing_step),
                      names_to = "stat",
                      values_to = "value"
                    ) %>%
                    # sort by processing_step factor levels
                    arrange(processing_step) %>%
                    # Create new column names by combining processing_step and stat
                    mutate(stat = paste(processing_step, stat, sep = "_")) %>%
                    dplyr::select(-processing_step) %>%
                    # Pivot back to wide format
                    pivot_wider(
                      names_from = stat,
                      values_from = value
                    ) 

  # If there is flT3_raw_total_sequences column, use it to calculate flT1_multimapping_reads
  if (any(grepl("flT3_raw_total_sequences", colnames(tmp)))) {
    actual_flt2_rts_colname <- paste0(opt$endogenous_genome_name, "_flT2_raw_total_sequences")
    actual_flt2_mmr_colname <- paste0(opt$endogenous_genome_name, "_flT2_multimapping_reads")
    actual_flt3_rts_colname <- paste0(opt$endogenous_genome_name, "_flT3_raw_total_sequences")

    # add new flT1_multimapping_reads or flT2_multimapping_reads column
    if (any(grepl(actual_flt2_rts_colname, colnames(tmp)))) {
      tmp <- tmp %>%
        mutate(!!sym(actual_flt2_mmr_colname) := !!sym(actual_flt2_rts_colname) - !!sym(actual_flt3_rts_colname)) %>%
        # move flT2_multimapping_reads after last column with flT2 in the name
        relocate(!!sym(actual_flt2_mmr_colname), .after = tail(grep("flT2", colnames(tmp[-which(colnames(tmp) == actual_flt2_mmr_colname)])), n = 1))
    } else if (any(grepl("flT1_raw_total_sequences", colnames(tmp)))) {
      tmp <- tmp %>%
        mutate(flT1_multimapping_reads = flT1_raw_total_sequences - !!sym(actual_flt3_rts_colname)) %>%
        # move flT1_multimapping_reads after last column with flT1 in the name
        relocate(flT1_multimapping_reads, .after = tail(grep("flT1", colnames(.[-which(colnames(.) == "flT1_multimapping_reads")])), n = 1))
    } else {
      warning("No valid columns found for flT1 or flT2 multimapping reads.")
      stop("Processing cannot continue without valid columns.")
    }
  }
}

# write the processed table to a file
write_tsv(tmp, file.path(opt$outdir, paste0(opt$prefix, ".tsv")))

cat("Processing completed on", Sys.time(), "\n")







                
                

                
