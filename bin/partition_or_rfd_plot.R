#!/usr/bin/env Rscript

# ===============================================================================
# partition_plot.R
#
# Originally created by:
#     - Nicolas Alcaraz <nicolas.alcaraz@cpr.ku.dk>
# Source:
#     https://github.com/grothlab/SCARseq_Pipeline/blob/a4b327f1901ae6a980767d05ec7af79896a604c9/SCAR_partition_plots.R
#
# Adapted for the grothlab/crepas pipeline by:
#     - Samuel Ruiz-Pérez <samper@cancer.dk>
#     https://github.com/grothlab/crepas/
#
# Description:
#     Script for plotting partition plots of SCARseq together with its stranded Input.
#     Also plots the RFD of one or more OK-seq files in the same figures and, for each
#     OK-seq file, scatter-correlation plots against the partitions.
# ===============================================================================

options(show.error.locations = TRUE)

required.libs <- c("tidyverse","GenomicAlignments","GenomicFeatures",
           "RColorBrewer","ggrepel","ggpubr","ggpmisc",
          "hexbin","argparse","sandwich","cowplot")

unavailable.libs <- setdiff(required.libs, rownames(installed.packages()))
if (length(unavailable.libs) > 0) {
  stop("Missing R packages: ", paste(unavailable.libs, collapse = ", "),
       ". Run the script with the conda environment or container of the PARTITION_OR_RFD_PLOT module.")
}

message("\n[", Sys.time(), "] Loading required libraries...")
for (lib in required.libs) {
  message("[", Sys.time(), "] Loading library: ", lib)
  suppressPackageStartupMessages({
    if (!require(lib, character.only = TRUE)) {
      stop("Failed to load library: ", lib)
    }
  })
}
message("[", Sys.time(), "] All libraries loaded successfully.")

# Plots built only to take parts of them (e.g. a legend) would otherwise be drawn to Rplots.pdf
grDevices::pdf(NULL)


# ===============================================================================
# Helper functions
# ===============================================================================

# Derive the sample name from a file name. Sample ids can themselves contain dots (e.g., H3.3).
# This function handles crepas extension conventions; --rfd_names overrides it for the RFD files.
sample_name_from_file <- function(path) {
  name <- basename(path)
  stripped <- sub("\\.(?:[0-9]+bp_[35]prime|cm_allo|allo|mLb|Lb)\\..*$", "", name, perl = TRUE)
  if (identical(stripped, name)) {
    stripped <- sub("\\.(bed|tsv|txt|bw|bigWig)(\\.gz)?$", "", name)
  }
  stripped
}


# Beyond 12 colors, the other qualitative Brewer palettes are appended
sample_palette <- function(n) {
  pool <- unique(c(
    RColorBrewer::brewer.pal(12, "Paired"),
    RColorBrewer::brewer.pal(8, "Dark2"),
    RColorBrewer::brewer.pal(9, "Set1"),
    RColorBrewer::brewer.pal(8, "Set2"),
    RColorBrewer::brewer.pal(12, "Set3")
  ))
  n_hues <- n
  while (length(pool) < n) {
    hues <- seq(15, 375, length.out = n_hues + 1)[seq_len(n_hues)]
    pool <- unique(c(pool,
                     grDevices::hcl(h = hues, c = 100,
                                    l = rep(c(65, 45, 80), length.out = n_hues))))
    n_hues <- n_hues * 2
  }
  pool[seq_len(n)]
}


# Line types that tell several RFD (OK-seq) lines apart, starting with a solid line
rfd_linetypes <- function(n) {
  pool <- c("solid", "22", "42", "13", "1343", "73", "2262",
            "12223242", "F282", "F4448444", "224282F2", "F1")
  if (n > length(pool)) {
    warning("\n[", Sys.time(), "] WARNING: More than ", length(pool),
            " RFD files were given; their line types are reused.")
  }
  rep_len(pool, n)
}


# Captions are drawn as the legend title, below the keys, one item per line
caption_lines <- function(...) {
  paste(c(...), collapse = "\n")
}


# Remove a leading "<exp_type>_" from sample names
strip_exp_type <- function(x, exp_type) {
  prefix <- paste0(exp_type, "_")
  ifelse(startsWith(x, prefix), substring(x, nchar(prefix) + 1), x)
}


# ===============================================================================
# Argument parsing
# ===============================================================================

parser <- ArgumentParser()

parser$add_argument("-a","--partition_files",
                    action = "store",
                    type = "character",
                    nargs = '+',
                    help = "Partition file(s), e.g. from SCAR-seq or eSPAN. It can be a single file path or a space-separated list of quoted file paths. If multiple files are provided, they will be plotted in the same figure.")

parser$add_argument("-d","--minusinput_partition_files", action = "store",
                    type = "character",
                    nargs = '+',
                    help = "Input-corrected partition file(s), in the same sample order as --partition_files. It can be a single file path or a space-separated list of quoted file paths. If multiple files are provided, they will be plotted in the same figure.")

parser$add_argument("-f","--strandedinput_partition_files", action = "store",
                    type = "character",
                    nargs = '+',
                    help = "Stranded-input partition file(s), in the same sample order as --partition_files. It can be a single file path or a space-separated list of quoted file paths. If multiple files are provided, they will be plotted in the same figure.")

parser$add_argument("-k", "--rfd_files", action = "store",
                    type = "character",
                    nargs = '+',
                    help = "RFD file(s) from OK-seq [optional, required for scatter plots]. It can be a single file path or a space-separated list of file paths. All of them are plotted together in the partition/RFD plots, with a different line type each, and each one gets its own set of scatter plots against the partitions.")

parser$add_argument("--rfd_exp_type", action = "store",
                    default = "OK-seq",
                    type = "character",
                    help = "Name of the experiment type of the RFD file(s), used in plot labels and legends. With several RFD files, each is labelled '<rfd_exp_type> (<sample name>)' [default: OK-seq]")

parser$add_argument("--rfd_names", action = "store",
                    type = "character",
                    help = "Comma-separated sample names of the RFD file(s), one per --rfd_files entry and in the same order [optional]. They replace the sample names derived from the file names in the plot legends, the scatter plot file names and the sample column of the mean values table. A leading '<rfd_exp_type>_' is removed from them, as from the derived names.")

parser$add_argument("--partition_exp_type", action = "store",
                    default = "SCAR-seq",
                    type = "character",
                    help = "Name of the experiment type of the partition files (e.g. SCAR-seq or eSPAN), used in the facet labels. A leading '<partition_exp_type>_' is removed from the sample names in the legends [default: SCAR-seq]")

parser$add_argument("-i", "--initiation_zones", action = "store",
                    type = "character",
                    required = TRUE,
                    help = "Bed file with known initiation-zones.")

parser$add_argument("-b", "--blacklist", action = "store",
                    type = "character",
                    help = "Blacklist bed file with regions to exclude [optional but strongly recommended]")

parser$add_argument("-s", "--chrom_sizes", action = "store",
                    type = "character",
                    help = "Chromosome sizes file [required]")

parser$add_argument("-n", "--prefix", action = "store",
                    default = "partition",
                    type = "character",
                    help = "Prefix for output files and plot title [default: partition]")

parser$add_argument("-o", "--outdir", action = "store",
                    default = "plots",
                    type = "character",
                    help = "Path to output directory for plots [default: ./plots]")

parser$add_argument("-c", "--rpm_cutoff", action = "store",
                    default = 0.3,
                    type = "double",
                    help = "RPM cutoff for noisy bins [default: 0.3]")

parser$add_argument("-r", "--plot_range", action = "store",
                    default = 100,
                    type = "integer",
                    help = "Distance (kb) up- and downstream of each initiation zone to plot. Independent of --iz_rm_overlap_range: when --plot_range is the larger of the two, the plotted flanks can reach neighbouring initiation zones and a partition bin may then contribute to more than one profile. [default: 100 kb]")

parser$add_argument("-z", "--iz_rm_overlap_range", action = "store",
                    default = 100,
                    type = "integer",
                    help = "Distance (kb) around each initiation zone used to discard overlapping ones. Any initiation zone whose window of this radius overlaps another one is discarded before plotting. Set to 0 to skip the removal of overlapping initiation zones entirely and keep them all. Independent of --plot_range, so the plotted range can be widened without discarding more initiation zones. [default: 100 kb]")

parser$add_argument("-e", "--exclude_chromosomes", action = "store",
                    default = "chrX,chrY,chrM",
                    type = "character",
                    help = "Chromosomes to exclude from analyses, must be provided comma separated [default: chrX,chrY,chrM]")

parser$add_argument("-g", "--exclude_scaffolds", action = "store",
                    default = TRUE,
                    type = "logical",
                    help = "Whether to exclude scaffolds from analyses. Chromosomes whose name begins with 'chrUn', ends with '_random', '_alt' or '_fix', or contains a dot ('.') are considered scaffolds [default: TRUE]")

parser$add_argument("-w", "--only_plot_wholly_within_iz", action = "store",
                    default = FALSE,
                    type = "logical",
                    help = "When calculating the distance from each partition/RFD bin to the start of initiation zones, only consider bins that are wholly contained within the initiation zones [default: FALSE]")

parser$add_argument("--violin_test", action = "store",
                    default = "wilcoxon",
                    choices = c("wilcoxon", "newey_west"),
                    type = "character",
                    help = "Paired test between samples in the violin plots, which compare one value per initiation zone: 'wilcoxon' (Wilcoxon signed-rank test) or 'newey_west' (t-test of the mean paired difference with a Newey-West standard error, with the initiation zones ordered along the genome, which allows for correlation between neighbouring initiation zones). P-values are Holm-adjusted [default: wilcoxon]")

parser$add_argument("--newey_west_lag", action = "store",
                    default = NULL,
                    type = "integer",
                    help = "Lag (in initiation zones) of the Newey-West standard error with --violin_test newey_west [default: chosen automatically, Newey & West (1994)]")

parser$add_argument("--profile_test_bin_size", action = "store",
                    default = 5,
                    type = "double",
                    help = "Bin size (kb) of the distance axis for the profile difference tests: the raw values of each initiation zone are averaged in bins of this size before comparing samples. Smaller bins are noisier and split broad differences into many small clusters, and bins smaller than the windows of the tables leave empty bins [default: 5 kb]")

parser$add_argument("--profile_test_permutations", action = "store",
                    default = 5000,
                    type = "integer",
                    help = "Number of sign-flip permutations of the cluster-based permutation test of the profile differences. The smallest possible p-value is 1 / (permutations + 1), multiplied by the number of pairs after adjustment [default: 5000]")

parser$add_argument("--profile_test_threshold", action = "store",
                    default = 0.05,
                    type = "double",
                    help = "Pointwise two-sided p-value that forms the clusters in the cluster-based permutation test: neighbouring distance bins whose paired t-test passes it (with the same sign) form a cluster [default: 0.05]")

parser$add_argument("--seed", action = "store",
                    default = 1,
                    type = "integer",
                    help = "Random seed of the permutation tests, so that their p-values are reproducible [default: 1]")

opt <- parser$parse_args()

opt_partition_files <- opt$partition_files
opt_minusinput_partition_files <- opt$minusinput_partition_files
opt_strandedinput_partition_files <- opt$strandedinput_partition_files
opt_rfd_files <- opt$rfd_files
opt_rfd_exp_type <- opt$rfd_exp_type
opt_rfd_names <- opt$rfd_names
opt_partition_exp_type <- opt$partition_exp_type
opt_initiation_zones <- opt$initiation_zones
opt_blacklist <- opt$blacklist
opt_chrom_sizes <- opt$chrom_sizes
opt_prefix <- opt$prefix
opt_outdir <- opt$outdir
opt_rpm_cutoff <- opt$rpm_cutoff
opt_plot_range <- opt$plot_range
opt_iz_rm_overlap_range <- opt$iz_rm_overlap_range
opt_exclude_chromosomes <- opt$exclude_chromosomes
opt_exclude_scaffolds <- opt$exclude_scaffolds
opt_only_plot_within_iz <- opt$only_plot_wholly_within_iz
opt_violin_test <- opt$violin_test
opt_newey_west_lag <- opt$newey_west_lag

opt_profile_test_bin_size <- opt$profile_test_bin_size
opt_profile_test_permutations <- opt$profile_test_permutations
opt_profile_test_threshold <- opt$profile_test_threshold
opt_seed <- opt$seed

if (!is.null(opt_newey_west_lag) && opt_newey_west_lag < 0) {
  stop("[", Sys.time(), "] ERROR: --newey_west_lag must be 0 or larger.")
}
if (opt_profile_test_bin_size <= 0) {
  stop("[", Sys.time(), "] ERROR: --profile_test_bin_size must be larger than 0.")
}
if (opt_profile_test_permutations < 1) {
  stop("[", Sys.time(), "] ERROR: --profile_test_permutations must be at least 1.")
}
if (opt_profile_test_threshold <= 0 || opt_profile_test_threshold >= 1) {
  stop("[", Sys.time(), "] ERROR: --profile_test_threshold must be between 0 and 1.")
}


# ===============================================================================
# Initialize partition files and flags
# ===============================================================================

part_files <- list("partition" = c(), "partition_input_corrected" = c(), "stranded_input_partition" = c())

# Check partition files
for (file in opt_partition_files) {
  if (!file.exists(file)) {
    stop("Partition file not found: ", file)
  }
  part_files[["partition"]] <- c(part_files[["partition"]], file)
}

# Check input-corrected partition files
for (file in opt_minusinput_partition_files) {
  if (!file.exists(file)) {
    stop("Input-corrected partition file not found: ", file)
  }
  part_files[["partition_input_corrected"]] <- c(part_files[["partition_input_corrected"]], file)
}

# Check stranded input partition files
for (file in opt_strandedinput_partition_files) {
  if (!file.exists(file)) {
    stop("Stranded input partition file not found: ", file)
  }
  part_files[["stranded_input_partition"]] <- c(part_files[["stranded_input_partition"]], file)
}

# Check OK-seq RFD files
rfd_files <- c()
rfd_found <- logical(0)
if (is.null(opt_rfd_files)) {
  warning("\n[", Sys.time(), "] RFD file not provided.")
} else {
  rfd_found <- file.exists(opt_rfd_files)
  for (file in opt_rfd_files[!rfd_found]) {
    warning("\n[", Sys.time(), "] RFD file not found: ", file)
  }
  rfd_files <- opt_rfd_files[rfd_found]
}
HAS_RFD <- length(rfd_files) > 0

# Check the RFD sample names: one per --rfd_files entry, none empty.
# The appended comma makes strsplit() keep a trailing empty name instead of dropping it.
if (!is.null(opt_rfd_names)) {
  rfd_given_names <- trimws(strsplit(paste0(opt_rfd_names, ","), ",", fixed = TRUE)[[1]])
  if (length(rfd_given_names) != length(opt_rfd_files)) {
    stop("[", Sys.time(), "] ERROR: --rfd_names has ", length(rfd_given_names), " name(s) but --rfd_files has ",
         length(opt_rfd_files), " file(s). Give one name per RFD file, in the same order.")
  }
  if (any(rfd_given_names == "")) {
    stop("[", Sys.time(), "] ERROR: --rfd_names contains an empty name: '", opt_rfd_names, "'")
  }
}

# Name each OK-seq file in the plots: a single one keeps the bare experiment type name,
# several are told apart by their sample names (from --rfd_names, or else from the file names).
# rfd_file_ids name their scatter plot files.
if (is.null(opt_rfd_names)) {
  rfd_base_names <- vapply(rfd_files, sample_name_from_file, character(1), USE.NAMES = FALSE)
} else {
  rfd_base_names <- rfd_given_names[rfd_found]
}
# Their names in the plots and in the scatter plot file names drop a leading "<rfd_exp_type>_"
rfd_short_names <- as.character(strip_exp_type(rfd_base_names, opt_rfd_exp_type))
if (length(rfd_files) == 0) {
  rfd_labels <- character(0)
} else if (length(rfd_files) == 1) {
  rfd_labels <- opt_rfd_exp_type
} else {
  rfd_labels <- make.unique(paste0(opt_rfd_exp_type, " (", rfd_short_names, ")"), sep = "_")
}
rfd_file_ids <- make.unique(gsub("[^A-Za-z0-9._-]", "_", rfd_short_names), sep = "_")

# Whether each OK-seq file has the same bin size as every partition file (needed for its scatter plots)
RFD_BIN_SIZES <- setNames(integer(length(rfd_labels)), rfd_labels)
RFD_BIN_SIZES_MATCH <- setNames(rep(TRUE, length(rfd_labels)), rfd_labels)

# Count how many types have at least one file
num_types_with_files <- sum(sapply(part_files, function(x) length(x) > 0))

if (num_types_with_files == 0 && !HAS_RFD) {
  stop("[", Sys.time(), "] ERROR: Please provide at least one partition file or RFD file to create plots.")
} else if (num_types_with_files == 0 && HAS_RFD) {
  # OK-seq only
  plot_width <- 6.2
  plot_suffix <- "RFD"
} else if (num_types_with_files == 1) {
  plot_width <- 6.2
  plot_suffix <- "partition"
} else if (num_types_with_files == 2) {
  plot_width <- 7.2
  plot_suffix <- "partition"
} else if (num_types_with_files == 3) {
  plot_width <- 8.2
  plot_suffix <- "partition"
}

# Check that all partition types with files have the same number of files
if (num_types_with_files > 1) {
  file_counts <- sapply(part_files, length)
  file_counts <- file_counts[file_counts > 0]

  if (length(unique(file_counts)) > 1) {
    stop("[", Sys.time(), "] ERROR: All partition types must have the same number of files. ",
         "Current counts: ", paste(names(file_counts), "=", file_counts, collapse = ", "))
  }
}

# Check initiation zones file (argparse only checks that --initiation_zones is given)
if (!file.exists(opt_initiation_zones)) {
  stop("[", Sys.time(), "] Initiation zones file not found: ", opt_initiation_zones)
}

# Check blacklist file
HAS_BLACKLIST <- FALSE
if (is.null(opt_blacklist)) {
  warning("\n[", Sys.time(), "] Blacklist file not provided, it is recommended to use a blacklist")
} else if (!file.exists(opt_blacklist)) {
  warning("\n[", Sys.time(), "] Blacklist file not found: ", opt_blacklist, ". It is recommended to use a blacklist.")
} else {
  HAS_BLACKLIST <- TRUE
}

# Check chromosome sizes file
if (is.null(opt_chrom_sizes)) {
  stop("[", Sys.time(), "] Chromosome sizes file not provided.")
} else if (!file.exists(opt_chrom_sizes)) {
  stop("[", Sys.time(), "] Chromosome sizes file not found: ", opt_chrom_sizes)
}

# Check output directory
if (!dir.exists(opt_outdir)) {
  warning("\n[", Sys.time(), "] Output directory not found, creating: ", opt_outdir)
  dir.create(opt_outdir, recursive = TRUE)
}

# Check other parameters
IZ_PLOT_LIMITS <- opt_plot_range * 1000
IZ_OVERLAP_LIMITS <- opt_iz_rm_overlap_range * 1000


# ===============================================================================
# Reading blacklist and chromosome sizes
# ===============================================================================

cls <- c("seqnames","start","end", # Coordinates of bin
         "fwd_counts","rev_counts",# Forward (fwd_counts) and Reverse (rev_counts) raw counts in bin
         "fwd_RPM","rev_RPM",      # Forward (fwd_RPM) and Reverse (rev_RPM) RPMs in bin
         "RFD_raw","RFD_smooth",   # Partition scores computed with raw (RFD_raw) and smoothed (RFD_smooth) RPMs
         "RFD_deriv",              # Value of the derivative of the partition at this bin
         "score",                  # not used
         "zero_deriv")             # second derivative at this bin

blacklist_gr <- GRanges()
if (HAS_BLACKLIST) {
  blacklist_df <- read_tsv(opt_blacklist, col_select = c(1:3), col_names = c("seqnames","start","end"), show_col_types = FALSE)
  blacklist_gr <- makeGRangesFromDataFrame(blacklist_df, starts.in.df.are.0based = TRUE)
  rm(blacklist_df)
}

chrom_sizes_df <- read_tsv(opt_chrom_sizes, col_select = c(1:2), col_names = c("chr", "sizes"), show_col_types = FALSE)

# Remove excluded chromosomes
chrom_excl <- unique(unlist(strsplit(opt_exclude_chromosomes, ",")))
message("\n[", Sys.time(), "] Excluding the following chromosomes: ", paste(chrom_excl, collapse = ", "))
chrom_sizes_df <- chrom_sizes_df[!chrom_sizes_df$chr %in% chrom_excl, ]

# Remove scaffolds from chrom_sizes if needed
if (opt_exclude_scaffolds) {
  message("\n[", Sys.time(), "] Removing scaffolds from chromosome sizes...")
  chrom_sizes_df <- chrom_sizes_df[!grepl("^chrUn|_random$|_alt$|_fix$|\\.", chrom_sizes_df$chr), ]
}

chrom_sizes <- deframe(chrom_sizes_df)

rm(chrom_sizes_df)


message("\n# ===============================================================================")
message("# STEP 1. Extracting initiation zones from provided BED file...")
message("# ===============================================================================")

iz_base_name <- sample_name_from_file(opt_initiation_zones)
message("\n[", Sys.time(), "] (", iz_base_name, ") Reading initiation zones file...")
IZ_df <- read_tsv(opt_initiation_zones,
                  col_select = c(1:3),
                  col_names = c("seqnames", "start", "end"),
                  show_col_types = FALSE)

message("\n[", Sys.time(), "] (", iz_base_name, ") Removing initiation zones within excluded chromosomes or outside of chrom_sizes...")
IZ_gr <- IZ_df %>%
  filter(seqnames %in% names(chrom_sizes)) %>%
  makeGRangesFromDataFrame(seqinfo = chrom_sizes,
                            keep.extra.columns = TRUE,
                            starts.in.df.are.0based = TRUE)

# We copy the interval now and not before with dplyr because the start
# coordinates are now 1-based thanks to starts.in.df.are.0based = TRUE
IZ_gr$interval <- paste0(seqnames(IZ_gr), ":", start(IZ_gr), "-", end(IZ_gr))

# Get original start and end coordinates for each initiation zone
IZ_gr$break_start <- start(IZ_gr)
IZ_gr$break_end <- end(IZ_gr)

# A radius of 0 switches overlap removal off, keeping every initiation zone
if (opt_iz_rm_overlap_range > 0) {

  message("\n[", Sys.time(), "] (", iz_base_name, ") Removing overlapping initiation zones (within ", opt_iz_rm_overlap_range, " kb upstream and ", opt_iz_rm_overlap_range, " kb downstream of another initiation zone)...")

  # Overlapping initiation zones are discarded based on their own radius, which is independent
  # of the plotted range: widening the plot does not discard more initiation zones.
  # Resizing can generate bins with negative start positions (out-of-bound), so we trim them
  iz_overlap_gr <- trim(resize(IZ_gr, IZ_OVERLAP_LIMITS * 2, fix = "center"))

  # Finding the nearest resized initiation zone to each resized initiation zone
  IZ_dist <- distanceToNearest(iz_overlap_gr)

  # Removing overlapping resized initiation zones
  overlapping_hits <- queryHits(subset(IZ_dist, IZ_dist@elementMetadata$distance == 0))
  # The following line removes overlapping IZs, and catches the case when there are no overlaps
  if (length(overlapping_hits) > 0) {
    IZ_gr <- IZ_gr[-overlapping_hits]
  } else {
    message("\n[", Sys.time(), "] (", iz_base_name, ") No overlapping initiation zones found within ", opt_iz_rm_overlap_range, " kb upstream and ", opt_iz_rm_overlap_range, " kb downstream of another initiation zone.")
  }

  # Remove temporary variables
  rm(iz_overlap_gr, IZ_dist, overlapping_hits)

} else {
  message("\n[", Sys.time(), "] (", iz_base_name, ") --iz_rm_overlap_range is 0, so overlapping initiation zones are not removed.")
}

# The retained initiation zones are then resized to the plotted range
IZ_gr <- trim(resize(IZ_gr, IZ_PLOT_LIMITS * 2, fix = "center"))

message("\n[", Sys.time(), "] (", iz_base_name, ") The number of initiation zones retained is: ", length(IZ_gr), ".")


# Fate of the windows of a table that lie within the plotted range of the initiation zones, for
# the window QC plot: used, below the RPM cutoff, blacklisted, or no reads. Tables often leave
# out the windows without reads, so those are counted as the expected windows (from the window
# size and step of the table) that are missing from it.
count_window_fates <- function(df, sample, facet) {
  df <- df %>%
    dplyr::filter(seqnames %in% names(chrom_sizes)) %>%
    dplyr::arrange(seqnames, start)
  gr <- makeGRangesFromDataFrame(df, seqinfo = chrom_sizes, keep.extra.columns = TRUE,
                                 starts.in.df.are.0based = TRUE)
  gr <- gr[if (opt_only_plot_within_iz) overlapsAny(gr, IZ_gr, type = "within") else overlapsAny(gr, IZ_gr)]
  has_reads <- gr$fwd_counts >= 1 | gr$rev_counts >= 1
  blacklisted <- has_reads & overlapsAny(gr, blacklist_gr, minoverlap = 1)
  below_cutoff <- has_reads & !blacklisted & (gr$fwd_RPM + gr$rev_RPM) < opt_rpm_cutoff

  # Expected windows: width w, starting every `step` bp from 0 on each chromosome
  widths <- table(df$end - df$start)
  w <- as.numeric(names(widths)[which.max(widths)])
  same_chrom <- df$seqnames[-1] == df$seqnames[-nrow(df)]
  steps <- table(diff(df$start)[same_chrom & diff(df$start) > 0])
  step <- if (length(steps) > 0) as.numeric(names(steps)[which.max(steps)]) else w
  ranges <- reduce(IZ_gr)
  a <- start(ranges)
  b <- end(ranges)
  chrom_len <- chrom_sizes[as.character(seqnames(ranges))]
  if (opt_only_plot_within_iz) {
    k_min <- ceiling((a - 1) / step)
    k_max <- floor((b - w) / step)
  } else {
    k_min <- pmax(0, ceiling((a - w) / step))
    k_max <- floor((b - 1) / step)
  }
  k_max <- pmin(k_max, floor((chrom_len - 1) / step))
  n_expected <- sum(pmax(0, k_max - k_min + 1))
  n_missing <- max(0, n_expected - length(gr))

  tibble(sample = sample, sample_facet = facet,
         fate = c("used", "below RPM cutoff", "blacklisted", "no reads"),
         n_windows = c(sum(has_reads & !blacklisted & !below_cutoff), sum(below_cutoff),
                       sum(blacklisted), sum(!has_reads) + n_missing))
}
window_qc_df <- tibble()


if (HAS_RFD) {

  message("\n# ===============================================================================")
  message("# STEP 2. Preprocessing of ", opt_rfd_exp_type, " (RFD) file(s)...")
  message("# ===============================================================================")

  RFD_df <- tibble()

  for (i in seq_along(rfd_files)) {

    rfd_file <- rfd_files[i]
    rfd_base_name <- rfd_base_names[i]
    rfd_label <- rfd_labels[i]

    message("\n[", Sys.time(), "] (", rfd_base_name, ") Reading ", opt_rfd_exp_type, " RFD file...")
    RFD_file_df <- read_tsv(rfd_file, col_names = cls, show_col_types = FALSE)

    if (ncol(RFD_file_df) != 12) {
      stop("[", Sys.time(), "] (", rfd_base_name, ") ERROR: The ", opt_rfd_exp_type, " RFD file is not in the correct format (ncol != 12)")
    }

    window_qc_df <- bind_rows(window_qc_df, count_window_fates(RFD_file_df, rfd_label, "RFD"))

    message("\n[", Sys.time(), "] (", rfd_base_name, ") Removing ", opt_rfd_exp_type, " bins within excluded chromosomes, outside of chrom_sizes, or with <1 fwd or <1 rev counts...")
    RFD_file_df <- RFD_file_df %>%
      dplyr::filter(seqnames %in% names(chrom_sizes),
                    (fwd_counts >= 1 | rev_counts >= 1)) %>%
      mutate(total_counts = fwd_counts + rev_counts, # total raw counts
             RPM = fwd_RPM + rev_RPM,                 # total counts per million
             sample = rfd_label,
             sample_type = opt_rfd_exp_type)

    message("\n[", Sys.time(), "] (", rfd_base_name, ") Converting the start positions of the ", opt_rfd_exp_type, " bins to 1-based...")
    RFD_file_gr <- makeGRangesFromDataFrame(RFD_file_df,
                                      seqinfo = chrom_sizes,
                                      keep.extra.columns = TRUE,
                                      starts.in.df.are.0based = TRUE)

    # We copy the interval now and not before with dplyr because the start
    # coordinates are now 1-based thanks to starts.in.df.are.0based = TRUE
    RFD_file_gr$interval <- paste0(seqnames(RFD_file_gr), ":", start(RFD_file_gr), "-", end(RFD_file_gr))

    message("\n[", Sys.time(), "] (", rfd_base_name, ") Calculating ", opt_rfd_exp_type, " bin size...")

    RFD_BIN_SIZES[rfd_label] <- width(RFD_file_gr)[1]

    message("\n[", Sys.time(), "] (", rfd_base_name, ") Removing ", opt_rfd_exp_type, " bins that overlap a blacklisted region...")
    if (HAS_BLACKLIST) {
      RFD_file_gr <- RFD_file_gr[!overlapsAny(RFD_file_gr, blacklist_gr, minoverlap = 1)]
    }

    # Finding out which initiation zones overlap which OK-seq bins
    if (opt_only_plot_within_iz) {
      overlap_pairs <- findOverlaps(query = RFD_file_gr, subject = IZ_gr, type = "within")
    } else {
      overlap_pairs <- findOverlaps(query = RFD_file_gr, subject = IZ_gr)
    }

    message("\n[", Sys.time(), "] (", rfd_base_name, ") Extracting the ", opt_rfd_exp_type, " bins that overlap with an initiation zone...")
    RFD_gr <- RFD_file_gr[queryHits(overlap_pairs)]

    # Adding to each OK-seq bin the interval of the initiation zone with which it overlaps
    RFD_gr$break_ID <- IZ_gr$interval[subjectHits(overlap_pairs)]

    message("\n[", Sys.time(), "] (", rfd_base_name, ") Calculating the distance from each ", opt_rfd_exp_type, " bin to the initiation zone with which it overlaps...")
    RFD_gr$dist <- start(RFD_gr) - IZ_gr$break_start[subjectHits(overlap_pairs)]

    RFD_df <- RFD_df %>%
      bind_rows(as_tibble(RFD_gr))

    # Remove temporary variables
    rm(RFD_file_df, RFD_file_gr, overlap_pairs, RFD_gr)

  }

} else {

  message("\n# ===============================================================================")
  message("# STEP 2. SKIPPED (No RFD file)")
  message("# ===============================================================================")

}

message("\n# ===============================================================================")
message("# STEP 3. Preprocessing of partition files...")
message("# ===============================================================================")

partition_df <- tibble()

for (type in names(part_files)) {

  message("\n[", Sys.time(), "] Processing partition files of type: ", type, "...")

  for (pair_id in seq_along(part_files[[type]])) {

    file <- part_files[[type]][pair_id]

    base_name <- sample_name_from_file(file)

    # A file given more than once (e.g. a stranded input shared by several samples) is read once:
    # its rows are copied for the scatter plot pairing, and left out of the profile means
    previous <- match(file, part_files[[type]][seq_len(pair_id - 1)])
    if (!is.na(previous)) {
      message("\n[", Sys.time(), "] (", base_name, ") Same file as number ", previous, " of type ", type, ", reusing it...")
      partition_df <- partition_df %>%
        bind_rows(partition_df %>%
                    dplyr::filter(sample_type == type, pair_id == previous) %>%
                    mutate(pair_id = .env$pair_id, duplicate = TRUE))
      next
    }

    message("\n[", Sys.time(), "] (", base_name, ") Reading partition file...")
    part_df <- read_tsv(file, col_names = cls, show_col_types = FALSE)

    if (ncol(part_df) != 12) {
      stop("[", Sys.time(), "] ERROR:", base_name, "is not a partition file format (ncol != 12)")
    }

    window_qc_df <- bind_rows(window_qc_df, count_window_fates(part_df, base_name, type))

    message("\n[", Sys.time(), "] (", base_name, ") Removing partition bins within excluded chromosomes, outside of chromosome sizes, or with <1 fwd or <1 rev counts...")
    part_df <- part_df %>%
      dplyr::filter(seqnames %in% names(chrom_sizes),
                  (fwd_counts >= 1 | rev_counts >= 1)) %>%
      mutate(total_counts = fwd_counts + rev_counts, # total raw counts
             RPM = fwd_RPM + rev_RPM,                # total counts pr. million
             sample = base_name,
             sample_type = type,
             sample_facet = type,
             pair_id = pair_id,
             duplicate = FALSE)

    part_gr <- makeGRangesFromDataFrame(part_df,
                                        seqinfo = chrom_sizes,
                                        keep.extra.columns = TRUE,
                                        starts.in.df.are.0based = TRUE)

    # We copy the interval now and not before with dplyr because the start
    # coordinates are now 1-based thanks to starts.in.df.are.0based = TRUE
    part_gr$interval <- paste0(seqnames(part_gr), ":", start(part_gr), "-", end(part_gr))

    message("\n[", Sys.time(), "] (", base_name, ") Removing partition bins that overlap a blacklisted region...")
    part_gr <- part_gr[!overlapsAny(part_gr, blacklist_gr, minoverlap = 1)]

    PART_BIN_SIZE <- width(part_gr)[1]
    message("\n[", Sys.time(), "] (", base_name, ") This partition's bin size is ", PART_BIN_SIZE, " bp.")

    for (rfd_label in names(RFD_BIN_SIZES)[RFD_BIN_SIZES != PART_BIN_SIZE]) {
      RFD_BIN_SIZES_MATCH[rfd_label] <- FALSE
      warning("\n[", Sys.time(), "] WARNING: Bin size of ", base_name,
              " (", PART_BIN_SIZE, " bp) is not the same bin size as in the provided RFD file of ", rfd_label, " (", RFD_BIN_SIZES[[rfd_label]], " bp). ",
              "The scatter plots against ", rfd_label, " pair bins by their coordinates and will be skipped; the profile plots are not affected.\n")
    }

  # As in Petryk et al. (2018; https://www-science.org/doi/10.1126/science.aau0294#supplementary-materials):
  message("\n[", Sys.time(), "] (", base_name, ") Calculating partition rates around initiation zones (", opt_plot_range, " kb upstream and ", opt_plot_range, " kb downstream of each IZ) by averaging values within each bin position...")

  # Finding out which initiation zones overlap which partition bins...
  if (opt_only_plot_within_iz) {
    overlap_pairs <- findOverlaps(query = part_gr, subject = IZ_gr, type = "within")
  } else {
    overlap_pairs <- findOverlaps(query = part_gr, subject = IZ_gr)
  }

  message("\n[", Sys.time(), "] (", base_name, ") Extracting the partition bins that overlap with an initiation zone...")
  partition_gr <- part_gr[queryHits(overlap_pairs)]

  # Adding to each partition bin the interval of the initiation zone with which it overlaps
  partition_gr$break_ID <- IZ_gr$interval[subjectHits(overlap_pairs)]

  message("\n[", Sys.time(), "] (", base_name, ") Calculating the distance from each partition bin to the initiation zone with which it overlaps...")
  partition_gr$dist <- start(partition_gr) - IZ_gr$break_start[subjectHits(overlap_pairs)]

  partition_df <- partition_df %>%
    bind_rows(as_tibble(partition_gr))

  message("\n[", Sys.time(), "] (", base_name, ") Finished processing partition file.")

  }

}

if (HAS_RFD) {

  message("\n[", Sys.time(), "] Expanding ", opt_rfd_exp_type, " RFD data for each partition (plotting purposes)...")
  # Expand RFD_df for each unique sample_facet in partition_mean_df
  sample_facets <- unique(partition_df$sample_facet)
  # if sample_facets is empty (there are no partitions, just RFD files):
  if (length(sample_facets) == 0 || is.null(sample_facets)) {
    sample_facets <- opt_rfd_exp_type
  }

  expanded_RFD <- tidyr::expand_grid(
    RFD_df,
    sample_facet = sample_facets
  )

  partition_df <- partition_df %>%
    bind_rows(expanded_RFD)

}

message("\n# ===============================================================================")
message("# STEP 4. Calculating mean partition/RFD rates around initiation zones...")
message("# ===============================================================================")
# As in Petryk et al. (2018; https://www-science.org/doi/10.1126/science.aau0294#supplementary-materials):

# Keep bins passing the RPM cutoff.
# Make sure copies of a file given more than once (duplicate = TRUE) are averaged only once
if (!"duplicate" %in% names(partition_df)) {
  partition_df$duplicate <- NA
}
partition_mean_df <- partition_df %>%
  dplyr::filter(RPM >= opt_rpm_cutoff,
                !(duplicate %in% TRUE))

# Get intersection of break_IDs between all samples (after RPM filtering)
common_break_IDs <- Reduce(intersect, split(partition_mean_df$break_ID, partition_mean_df$sample))

partition_mean_df <- partition_mean_df %>%
  dplyr::group_by(dist, sample, sample_type, sample_facet) %>%
  dplyr::summarise(
    RFD_raw_mean = mean(RFD_raw, na.rm = TRUE),
    RFD_raw_sd = sd(RFD_raw, na.rm = TRUE),
    RFD_raw_n = sum(!is.na(RFD_raw)),
    RFD_smooth_mean = mean(RFD_smooth, na.rm = TRUE),
    RFD_smooth_sd = sd(RFD_smooth, na.rm = TRUE),
    RFD_smooth_n = sum(!is.na(RFD_smooth)),
    .groups = "drop"
  ) %>%
  mutate(sample = ifelse(sample %in% rfd_labels, sample, strip_exp_type(sample, opt_partition_exp_type)))

write_tsv(partition_mean_df,
           file = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_mean_values.tsv")),
           col_names = TRUE)

message("\n[", Sys.time(), "] A glimpse of the partition mean data frame:")
print(partition_mean_df)

# Medians and quartiles across initiation zones, for the interquartile range plots
partition_quantile_df <- partition_df %>%
  dplyr::filter(RPM >= opt_rpm_cutoff,
                !(duplicate %in% TRUE)) %>%
  dplyr::group_by(dist, sample, sample_type, sample_facet) %>%
  dplyr::summarise(
    RFD_raw_median = median(RFD_raw, na.rm = TRUE),
    RFD_raw_q25 = quantile(RFD_raw, 0.25, na.rm = TRUE, names = FALSE),
    RFD_raw_q75 = quantile(RFD_raw, 0.75, na.rm = TRUE, names = FALSE),
    RFD_raw_n = sum(!is.na(RFD_raw)),
    RFD_smooth_median = median(RFD_smooth, na.rm = TRUE),
    RFD_smooth_q25 = quantile(RFD_smooth, 0.25, na.rm = TRUE, names = FALSE),
    RFD_smooth_q75 = quantile(RFD_smooth, 0.75, na.rm = TRUE, names = FALSE),
    RFD_smooth_n = sum(!is.na(RFD_smooth)),
    .groups = "drop"
  ) %>%
  mutate(sample = ifelse(sample %in% rfd_labels, sample, strip_exp_type(sample, opt_partition_exp_type)))

write_tsv(partition_quantile_df,
          file = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_quantile_values.tsv")),
          col_names = TRUE)

message("\n# ===============================================================================")
message("# STEP 5. Plotting")
message("# ===============================================================================")


# ===============================================================================
# Partition plots labels and colors
# ===============================================================================

message("\n[", Sys.time(), "] Setting sample labels and colors...")

# Check if there are partition files
if (num_types_with_files > 0) {
  sample_labels <- c("partition_input_corrected" = paste0(opt_partition_exp_type, " (Input-corrected)"),
                     "partition" = opt_partition_exp_type, "stranded_input_partition" = "Stranded input")
  sample_labels <- sample_labels[names(part_files)]
} else if (HAS_RFD) {
  # Only OK-seq, no partition files
  sample_labels <- setNames(opt_rfd_exp_type, opt_rfd_exp_type)
}

# set a color for each sample, but set OK-seq to dark gray
sample_names <- unique(partition_mean_df$sample)
sample_colors <- setNames(sample_palette(length(sample_names)), sample_names)

# A legend with many samples takes up more of the figure, so widen it to keep the panel readable
if (length(sample_names) > 15) {
  plot_width <- plot_width + 3
}
rfd_plotted <- intersect(rfd_labels, sample_names)
# Only set OK-seq to grey if there are partition files
if (num_types_with_files > 0) {
  sample_colors[rfd_plotted] <- "grey60"
}
line_colors <- sample_colors

# Partition lines are solid; several OK-seq lines get a line type each
line_types <- setNames(rep("solid", length(sample_names)), sample_names)
if (length(rfd_plotted) > 1) {
  line_types[rfd_plotted] <- rfd_linetypes(length(rfd_plotted))
}
# Dashed legend keys need more width to be told apart
legend_key_width <- unit(ifelse(length(rfd_plotted) > 1, 0.8, 0.3), "cm")


# ===============================================================================
# Smoothed partition plot(s)
# ===============================================================================

message("\n[", Sys.time(), "] Creating smoothed partition plot(s)...")


# band: "none" (mean), "sd" (mean +/- 1 SD) or "iqr" (median and interquartile range)
plot_partition <- function(value, band = c("none", "sd", "iqr")) {

  band <- match.arg(band)
  if (band == "iqr") {
    plot_df <- partition_quantile_df
    center_col <- paste0(value, "_median")
    caption <- caption_lines(paste0("N = ", length(IZ_gr), " initiation zones"),
                             "Line: median",
                             "Band: interquartile range (25th-75th percentile)")
  } else {
    plot_df <- partition_mean_df
    center_col <- paste0(value, "_mean")
    caption <- if (band == "sd") {
      caption_lines(paste0("N = ", length(IZ_gr), " initiation zones"),
                    "Line: mean",
                    "Ribbon: +/- 1 SD across initiation zones")
    } else {
      caption_lines(paste0("N = ", length(IZ_gr), " initiation zones"),
                    "Line: mean")
    }
  }

  p <- ggplot(plot_df, aes(x = dist / 1000, y = .data[[center_col]], color = sample, linetype = sample)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = "grey96") +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             fill = "grey96") +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3)

  if (band == "sd") {
    p <- p +
      geom_ribbon(aes(ymin = .data[[center_col]] - .data[[paste0(value, "_sd")]],
                      ymax = .data[[center_col]] + .data[[paste0(value, "_sd")]], fill = sample),
                  linetype = 0,
                  alpha = 0.2)
  } else if (band == "iqr") {
    p <- p +
      geom_ribbon(aes(ymin = .data[[paste0(value, "_q25")]],
                      ymax = .data[[paste0(value, "_q75")]], fill = sample),
                  linetype = 0,
                  alpha = 0.2)
  }

  p +
    geom_line(linewidth = 0.3) +
    scale_color_manual(values = line_colors) +
    scale_linetype_manual(values = line_types) +
    scale_fill_manual(values = line_colors) +
    xlab("Distance from initiation zone center (kb)") +
    ylab(ifelse(HAS_RFD, ifelse(num_types_with_files > 0, "Partition or RFD", "RFD"), "Partition")) +
    do.call(labs, c(list(color = caption, linetype = caption), if (band != "none") list(fill = caption))) +
    theme_bw(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 5),
          legend.title.position = "bottom",
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          legend.key.width = legend_key_width,
          legend.key = element_rect(color = NA, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(size = 7),
          aspect.ratio = 1) +
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1),
           fill = guide_legend(order = 1)) +
    facet_wrap(~ sample_facet, nrow = 1, labeller = labeller(.cols = sample_labels)) +
    scale_y_continuous(expand = expansion(mult = c(0.11, 0.11))) +
    annotate(geom = 'text', label = 'Lagging', x = -Inf, y = Inf, hjust = -0.13,
            vjust = 2, size = 2) +
    annotate(geom = 'text', label = 'Leading', x = -Inf, y = -Inf, hjust = -0.13,
            vjust = -1.5, size = 2) +
    annotate(geom = 'text', label = 'Lagging', x = Inf, y = -Inf, hjust = 1.13,
            vjust = -1.5, size = 2) +
    annotate(geom = 'text', label = 'Leading', x = Inf, y = Inf, hjust = 1.13,
            vjust = 2, size = 2)
}

smoothed_plot     <- plot_partition("RFD_smooth", "none")
smoothed_plot_sd  <- plot_partition("RFD_smooth", "sd")
smoothed_plot_iqr <- plot_partition("RFD_smooth", "iqr")

message("\n[", Sys.time(), "] Saving smoothed partition plot(s)...")
ggsave(filename = file.path(opt_outdir,paste0(opt_prefix,".", plot_suffix, "_plot_smoothed.pdf")),
       plot = smoothed_plot, width = plot_width, height = 2.5, units = "in")

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_smoothed.png")),
       plot = smoothed_plot, width = plot_width, height = 2.5, units = "in",
       dpi = 600)

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_smoothed.sd.pdf")),
       plot = smoothed_plot_sd, width = plot_width, height = 2.5, units = "in")

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_smoothed.sd.png")),
       plot = smoothed_plot_sd, width = plot_width, height = 2.5, units = "in",
       dpi = 600)

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_smoothed.iqr.pdf")),
       plot = smoothed_plot_iqr, width = plot_width, height = 2.5, units = "in")

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_smoothed.iqr.png")),
       plot = smoothed_plot_iqr, width = plot_width, height = 2.5, units = "in",
       dpi = 600)


# ===============================================================================
# Raw (unsmoothed) partition plot(s)
# ===============================================================================

message("\n[", Sys.time(), "] Creating raw partition plot(s)...")

raw_plot     <- plot_partition("RFD_raw", "none")
raw_plot_sd  <- plot_partition("RFD_raw", "sd")
raw_plot_iqr <- plot_partition("RFD_raw", "iqr")

message("\n[", Sys.time(), "] Saving raw partition plot(s)...")
ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_raw.pdf")),
       plot = raw_plot, width = plot_width, height = 2.5, units = "in")

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_raw.png")),
       plot = raw_plot, width = plot_width, height = 2.5, units = "in",
       dpi = 600)

message("\n[", Sys.time(), "] Saving raw partition plot(s) with the standard deviation ribbon...")
ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_raw.sd.pdf")),
       plot = raw_plot_sd, width = plot_width, height = 2.5, units = "in")

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_raw.sd.png")),
       plot = raw_plot_sd, width = plot_width, height = 2.5, units = "in",
       dpi = 600)

message("\n[", Sys.time(), "] Saving raw partition plot(s) with the interquartile band...")
ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_raw.iqr.pdf")),
       plot = raw_plot_iqr, width = plot_width, height = 2.5, units = "in")

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_raw.iqr.png")),
       plot = raw_plot_iqr, width = plot_width, height = 2.5, units = "in",
       dpi = 600)


# ===============================================================================
# GAM-fitted smoothed partition plot(s)
# ===============================================================================

message("\n[", Sys.time(), "] Creating GAM-fitted partition plot(s)...")

smoothed_gam_plot <- ggplot(partition_mean_df, aes(x = dist / 1000, y = RFD_smooth_mean, color = sample, linetype = sample)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = "grey96") +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             fill = "grey96") +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(stat = "smooth", method = "gam", se = FALSE, linewidth = 0.5) +
    scale_color_manual(values = line_colors) +
    scale_linetype_manual(values = line_types) +
    xlab("Distance from initiation zone center (kb)") +
    ylab(ifelse(HAS_RFD, ifelse(num_types_with_files > 0, "Partition or RFD", "RFD"), "Partition")) +
    labs(color = caption_lines(paste0("N = ", length(IZ_gr), " initiation zones"), "Line: GAM fit of the mean"),
         linetype = caption_lines(paste0("N = ", length(IZ_gr), " initiation zones"), "Line: GAM fit of the mean")) +
    theme_bw(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 5),
          legend.title.position = "bottom",
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          legend.key.width = legend_key_width,
          legend.key = element_rect(color = NA, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(size = 7),
          aspect.ratio = 1) +
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1)) +
    facet_wrap(~ sample_facet, nrow = 1, labeller = labeller(.cols = sample_labels)) +
    scale_y_continuous(expand = expansion(mult = c(0.11, 0.11))) +
    annotate(geom = 'text', label = 'Lagging', x = -Inf, y = Inf, hjust = -0.13,
            vjust = 2, size = 2) +
    annotate(geom = 'text', label = 'Leading', x = -Inf, y = -Inf, hjust = -0.13,
            vjust = -1.5, size = 2) +
    annotate(geom = 'text', label = 'Lagging', x = Inf, y = -Inf, hjust = 1.13,
            vjust = -1.5, size = 2) +
    annotate(geom = 'text', label = 'Leading', x = Inf, y = Inf, hjust = 1.13,
            vjust = 2, size = 2)

message("\n[", Sys.time(), "] Saving smoothed partition plot(s)...")
ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_smoothed.gam.pdf")),
       plot = smoothed_gam_plot, width = plot_width, height = 2.5, units = "in")

ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_plot_smoothed.gam.png")),
       plot = smoothed_gam_plot, width = plot_width, height = 2.5, units = "in",
       dpi = 600)


# ===============================================================================
# Violin plots: one value per initiation zone and sample, with paired tests between samples
# ===============================================================================
# Values change sign across an initiation zone, so each one gets
# (downstream mean - upstream mean) / 2 of its raw values; bins inside the zone are left out.

message("\n[", Sys.time(), "] Creating violin plot(s) of the per-initiation-zone values...")

# Wilcoxon signed-rank test of paired values; NA when it cannot be computed
wilcoxon_paired_p <- function(x, y) {
  if (length(x) < 2 || all(x - y == 0)) return(NA_real_)
  tryCatch(suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value),
           error = function(e) NA_real_)
}

# Test of mean(d) = 0 with a Newey-West standard error (d ordered along the genome)
newey_west_p <- function(d, lag) {
  if (length(d) < 3 || all(d == d[1])) return(NA_real_)
  fit <- lm(d ~ 1)
  se <- sqrt(sandwich::NeweyWest(fit, lag = lag, prewhite = FALSE)[1, 1])
  if (!is.finite(se) || se == 0) return(NA_real_)
  2 * pnorm(-abs(unname(coef(fit)[1]) / se))
}

p_to_stars <- function(p) {
  dplyr::case_when(is.na(p)     ~ "NA",
                   p <= 0.0001  ~ "****",
                   p <= 0.001   ~ "***",
                   p <= 0.01    ~ "**",
                   p <= 0.05    ~ "*",
                   TRUE         ~ "ns")
}

# Original (not resized) initiation zone coordinates, to tell the bins on either side apart
iz_info <- tibble(break_ID = IZ_gr$interval,
                  iz_seqnames = as.character(seqnames(IZ_gr)),
                  iz_start = IZ_gr$break_start,
                  iz_end = IZ_gr$break_end)

display_name <- function(x) ifelse(x %in% rfd_labels, x, strip_exp_type(x, opt_partition_exp_type))

# Samples in input order: partitions first (by type, then file order), then the RFD tables
violin_sample_levels <- unique(display_name(unique(partition_df$sample)))

iz_values_df <- partition_df %>%
  dplyr::filter(RPM >= opt_rpm_cutoff,
                !(duplicate %in% TRUE),
                !is.na(RFD_raw)) %>%
  dplyr::inner_join(iz_info, by = "break_ID") %>%
  mutate(side = dplyr::case_when(end < iz_start ~ "upstream",
                                 start > iz_end ~ "downstream",
                                 TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(side)) %>%
  group_by(sample, sample_type, sample_facet, break_ID, iz_seqnames, iz_start) %>%
  summarise(n_upstream = sum(side == "upstream"),
            n_downstream = sum(side == "downstream"),
            value = (mean(RFD_raw[side == "downstream"]) - mean(RFD_raw[side == "upstream"])) / 2,
            .groups = "drop") %>%
  dplyr::filter(n_upstream > 0, n_downstream > 0) %>%
  mutate(sample = display_name(sample),
         is_rfd = sample %in% rfd_labels) %>%
  dplyr::arrange(sample_facet, factor(sample, levels = violin_sample_levels), iz_seqnames, iz_start)

# With partitions, the RFD tables (repeated in every partition facet) get a single facet of their own
violin_facet_labels <- sample_labels
if (num_types_with_files > 0 && any(iz_values_df$is_rfd)) {
  first_rfd_facet <- unique(iz_values_df$sample_facet[iz_values_df$is_rfd])[1]
  iz_values_df <- iz_values_df %>%
    dplyr::filter(!is_rfd | sample_facet == first_rfd_facet) %>%
    mutate(sample_facet = ifelse(is_rfd, "RFD", sample_facet))
  violin_facet_labels <- c(sample_labels, RFD = opt_rfd_exp_type)
}

write_tsv(iz_values_df %>% dplyr::select(sample, sample_type, sample_facet, break_ID,
                                         n_upstream, n_downstream, value),
          file = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_violin_values.tsv")),
          col_names = TRUE)

# Paired tests between all pairs of partitions and all pairs of RFD tables in each facet
violin_stats_df <- tibble(sample_facet = character(), sample_kind = character(),
                          group1 = character(), group2 = character(), n_iz = integer(),
                          median1 = double(), median2 = double(), median_difference = double(),
                          mean_difference = double(), test = character(),
                          p_value = double(), p_adj = double(), p_adj_signif = character())

for (facet in unique(iz_values_df$sample_facet)) {
  for (kind in c(FALSE, TRUE)) {
    facet_df <- iz_values_df %>% dplyr::filter(sample_facet == facet, is_rfd == kind)
    samples <- intersect(violin_sample_levels, facet_df$sample)
    if (length(samples) < 2) next
    for (pair in combn(samples, 2, simplify = FALSE)) {
      paired <- dplyr::inner_join(
        facet_df %>% dplyr::filter(sample == pair[1]) %>% dplyr::select(break_ID, iz_seqnames, iz_start, x = value),
        facet_df %>% dplyr::filter(sample == pair[2]) %>% dplyr::select(break_ID, y = value),
        by = "break_ID") %>%
        dplyr::arrange(iz_seqnames, iz_start)
      d <- paired$x - paired$y
      p <- if (opt_violin_test == "wilcoxon") {
        wilcoxon_paired_p(paired$x, paired$y)
      } else {
        newey_west_p(d, opt_newey_west_lag)
      }
      violin_stats_df <- violin_stats_df %>%
        bind_rows(tibble(sample_facet = facet,
                         sample_kind = ifelse(kind, "RFD", "partition"),
                         group1 = pair[1], group2 = pair[2], n_iz = nrow(paired),
                         median1 = median(paired$x), median2 = median(paired$y),
                         median_difference = median(d), mean_difference = mean(d),
                         test = opt_violin_test, p_value = p))
    }
  }
}

violin_stats_df <- violin_stats_df %>%
  group_by(sample_facet, sample_kind) %>%
  mutate(p_adj = p.adjust(p_value, method = "holm")) %>%
  ungroup() %>%
  mutate(p_adj_signif = p_to_stars(p_adj))

write_tsv(violin_stats_df,
          file = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_violin_stats.tsv")),
          col_names = TRUE, na = "NA")

if (nrow(iz_values_df) == 0) {

  message("\n[", Sys.time(), "] Skipping the violin plot(s): no initiation zone has bins on both sides in any sample.")

} else {

  # Partition facets in their usual order, then the RFD facet
  violin_facet_levels <- c(sort(setdiff(unique(iz_values_df$sample_facet), "RFD")),
                           intersect("RFD", unique(iz_values_df$sample_facet)))
  iz_values_df <- iz_values_df %>%
    mutate(sample = factor(sample, levels = violin_sample_levels),
           sample_facet = factor(sample_facet, levels = violin_facet_levels))

  # Brackets: x positions within each facet (free x scales keep only that facet's samples),
  # stacked in tiers so that overlapping brackets do not cross; heights are set per plot below
  brackets_df <- tibble()
  for (facet in unique(violin_stats_df$sample_facet)) {
    facet_levels <- intersect(violin_sample_levels, as.character(iz_values_df$sample[iz_values_df$sample_facet == facet]))
    facet_stats <- violin_stats_df %>%
      dplyr::filter(sample_facet == facet) %>%
      mutate(xmin = match(group1, facet_levels), xmax = match(group2, facet_levels)) %>%
      dplyr::arrange(xmax - xmin, xmin)
    tiers <- integer(nrow(facet_stats))
    occupied <- list()
    for (i in seq_len(nrow(facet_stats))) {
      tier <- 1
      while (any(vapply(occupied[[as.character(tier)]] %||% list(),
                        function(span) facet_stats$xmin[i] <= span[2] && facet_stats$xmax[i] >= span[1],
                        logical(1)))) {
        tier <- tier + 1
      }
      occupied[[as.character(tier)]] <- c(occupied[[as.character(tier)]] %||% list(),
                                          list(c(facet_stats$xmin[i], facet_stats$xmax[i])))
      tiers[i] <- tier
    }
    brackets_df <- bind_rows(brackets_df, facet_stats %>% mutate(tier = tiers))
  }
  n_tiers <- if (nrow(brackets_df) > 0) max(brackets_df$tier) else 0

  # Median line in the violins; ggplot2 4.0 replaced draw_quantiles by quantiles
  # The RFD tables share one grey; as in the profile plots, their outlines get a line type each
  violin_layer <- if (utils::packageVersion("ggplot2") >= "4.0.0") {
    geom_violin(aes(fill = sample, linetype = sample), quantiles = 0.5, quantile.linetype = "solid",
                linewidth = 0.3, alpha = 0.6, scale = "width", trim = TRUE)
  } else {
    geom_violin(aes(fill = sample, linetype = sample), draw_quantiles = 0.5,
                linewidth = 0.3, alpha = 0.6, scale = "width", trim = TRUE)
  }

  n_iz_used <- n_distinct(iz_values_df$break_ID)
  test_name <- ifelse(opt_violin_test == "wilcoxon", "Wilcoxon signed-rank test",
                      paste0("Newey-West test (lag ", ifelse(is.null(opt_newey_west_lag), "automatic", opt_newey_west_lag), ")"))

  violin_caption <- caption_lines(paste0("N = ", n_iz_used, " initiation zones"),
                                  paste0("Test: ", test_name, ", paired by initiation zone"),
                                  "P-values: Holm-adjusted over the pairs")
  max_label_chars <- max(nchar(as.character(iz_values_df$sample)))
  n_per_facet <- iz_values_df %>% dplyr::distinct(sample_facet, sample) %>% dplyr::count(sample_facet)

  # One y scale for all the violins, with room for the brackets
  y_max <- max(iz_values_df$value, na.rm = TRUE)
  y_range <- diff(range(iz_values_df$value, na.rm = TRUE))
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  tip <- 0.02 * y_range
  violin_ylim <- c(min(iz_values_df$value, na.rm = TRUE),
                   y_max + (n_tiers + 0.9) * 0.09 * y_range)

  # Violin plot of the given facets, with its own y axis title
  violin_subplot <- function(facets, value_name) {
    df <- iz_values_df %>%
      dplyr::filter(sample_facet %in% facets) %>%
      mutate(sample_facet = factor(as.character(sample_facet), levels = intersect(violin_facet_levels, facets)))
    label_chars <- max(nchar(as.character(df$sample)))

    p <- ggplot(df, aes(x = sample, y = value)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
      violin_layer +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
      scale_fill_manual(values = sample_colors) +
      scale_linetype_manual(values = line_types) +
      xlab(NULL) +
      ylab(paste0("(mean downstream raw ", value_name, " -\nmean upstream raw ", value_name,
                  ") / 2, +/- ", opt_plot_range, " kb")) +
      labs(fill = violin_caption, linetype = violin_caption) +
      theme_bw(base_family = "Helvetica") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_text(size = 5),
            legend.title.position = "bottom",
            legend.text = element_text(size = 6),
            legend.key.size = unit(0.3, "cm"),
            legend.key.width = legend_key_width,
            axis.title = element_text(size = 6),
            axis.text = element_text(size = 6),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_blank(),
            strip.text = element_text(size = 7),
            # The first x label is rotated to the left of the panel: widen the margin for long labels
            plot.margin = margin(5.5, 5.5, 5.5, max(5.5, 2.2 * label_chars - 35), unit = "pt")) +
      facet_wrap(~ sample_facet, nrow = 1, scales = "free_x", labeller = labeller(.cols = violin_facet_labels)) +
      expand_limits(y = violin_ylim)

    br <- if (nrow(brackets_df) > 0) brackets_df %>% dplyr::filter(sample_facet %in% facets) else brackets_df
    if (nrow(br) > 0) {
      br <- br %>%
        mutate(y = y_max + tier * 0.09 * y_range,
               sample_facet = factor(as.character(sample_facet), levels = levels(df$sample_facet)))
      p <- p +
        geom_segment(data = br, aes(x = xmin, xend = xmax, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = 0.3) +
        geom_segment(data = br, aes(x = xmin, xend = xmin, y = y, yend = y - tip),
                     inherit.aes = FALSE, linewidth = 0.3) +
        geom_segment(data = br, aes(x = xmax, xend = xmax, y = y, yend = y - tip),
                     inherit.aes = FALSE, linewidth = 0.3) +
        geom_text(data = br, aes(x = (xmin + xmax) / 2, y = y, label = p_adj_signif),
                  inherit.aes = FALSE, vjust = -0.3, size = 2)
    }
    p
  }

  # Width of the panels of the given facets, with the y axis and the margin for rotated labels
  subplot_width <- function(facets) {
    label_chars <- max(nchar(as.character(iz_values_df$sample[iz_values_df$sample_facet %in% facets])))
    0.9 + sum(0.6 + 0.45 * n_per_facet$n[n_per_facet$sample_facet %in% facets]) + max(0, 0.031 * label_chars - 0.49)
  }
  legend_width <- 0.9 + max(0.05 * max_label_chars, 1.7)

  # With partitions and RFD tables, each gets its own plot (and y axis title), side by side
  partition_facets <- setdiff(violin_facet_levels, "RFD")
  if (length(partition_facets) > 0 && "RFD" %in% violin_facet_levels) {
    widths <- c(subplot_width(partition_facets), subplot_width("RFD"), legend_width)
    violin_legend <- ggpubr::get_legend(violin_subplot(violin_facet_levels, "partition"))
    # Panels aligned at top and bottom, whatever the height of their sample labels
    violin_plot <- cowplot::plot_grid(violin_subplot(partition_facets, "partition") + theme(legend.position = "none"),
                                      violin_subplot("RFD", "RFD") + theme(legend.position = "none"),
                                      violin_legend,
                                      nrow = 1, rel_widths = widths, align = "h", axis = "tb")
  } else {
    widths <- c(subplot_width(violin_facet_levels), legend_width)
    violin_plot <- violin_subplot(violin_facet_levels, ifelse(length(partition_facets) > 0, "partition", "RFD"))
  }
  violin_width <- sum(widths)
  # The 45-degree sample names take about 0.03 in of height per character
  violin_height <- 2.8 + 0.18 * n_tiers + 0.03 * max_label_chars

  message("\n[", Sys.time(), "] Saving violin plot(s)...")
  ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_violin_plot_raw.pdf")),
         plot = violin_plot, width = violin_width, height = violin_height, units = "in",
         limitsize = FALSE)

  ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_violin_plot_raw.png")),
         plot = violin_plot, width = violin_width, height = violin_height, units = "in",
         dpi = 600, limitsize = FALSE, bg = "white")

}


# ===============================================================================
# Window QC: fate of the windows around the initiation zones, and coverage along the distance axis
# ===============================================================================

message("\n[", Sys.time(), "] Creating the window QC plot...")

# As in the violin plots, the RFD tables get a facet of their own when there are partitions
qc_rfd_facet <- if (num_types_with_files > 0) "RFD" else opt_rfd_exp_type
qc_facet_labels <- c(sample_labels, RFD = opt_rfd_exp_type)

window_qc_df <- window_qc_df %>%
  mutate(sample = display_name(sample),
         sample_facet = ifelse(sample_facet == "RFD", qc_rfd_facet, sample_facet)) %>%
  group_by(sample, sample_facet) %>%
  mutate(percent = 100 * n_windows / sum(n_windows)) %>%
  ungroup()

write_tsv(window_qc_df,
          file = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_window_qc.tsv")),
          col_names = TRUE)

# Percentage of initiation zones with a used window at each distance (RFD tables counted once)
coverage_df <- partition_df %>%
  dplyr::filter(RPM >= opt_rpm_cutoff,
                !(duplicate %in% TRUE),
                !is.na(RFD_raw)) %>%
  mutate(sample = display_name(sample),
         is_rfd = sample %in% rfd_labels)
if (any(coverage_df$is_rfd)) {
  first_rfd_facet <- unique(coverage_df$sample_facet[coverage_df$is_rfd])[1]
  coverage_df <- coverage_df %>%
    dplyr::filter(!is_rfd | sample_facet == first_rfd_facet) %>%
    mutate(sample_facet = ifelse(is_rfd, qc_rfd_facet, sample_facet))
}
coverage_df <- coverage_df %>%
  group_by(sample, sample_facet, dist) %>%
  summarise(percent_iz = 100 * n_distinct(break_ID) / length(IZ_gr), .groups = "drop")

if (nrow(window_qc_df) > 0 && sum(window_qc_df$n_windows) > 0) {

  qc_facet_levels <- c(sort(setdiff(unique(window_qc_df$sample_facet), qc_rfd_facet)),
                       intersect(qc_rfd_facet, unique(window_qc_df$sample_facet)))
  qc_sample_levels <- unique(c(intersect(violin_sample_levels, window_qc_df$sample), window_qc_df$sample))
  fate_colors <- c("used" = "#4D4D4D", "below RPM cutoff" = "#E8C547",
                   "blacklisted" = "#C0504D", "no reads" = "#BFBFBF")

  qc_theme <- theme_bw(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 5),
          legend.title.position = "bottom",
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          strip.background = element_blank(),
          strip.text = element_text(size = 7))

  # A: fate of the windows within the plotted range, per sample
  qc_fates_plot <- ggplot(window_qc_df %>%
                            mutate(sample = factor(sample, levels = qc_sample_levels),
                                   sample_facet = factor(sample_facet, levels = qc_facet_levels),
                                   fate = factor(fate, levels = rev(names(fate_colors)))),
                          aes(x = sample, y = percent, fill = fate)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = fate_colors, guide = guide_legend(reverse = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    facet_wrap(~ sample_facet, nrow = 1, scales = "free_x", labeller = labeller(.cols = qc_facet_labels)) +
    xlab(NULL) +
    ylab(paste0("Windows within +/- ", opt_plot_range, " kb (%)")) +
    labs(fill = caption_lines(paste0("Windows within +/- ", opt_plot_range, " kb of ",
                                     length(IZ_gr), " initiation zones"),
                              paste0("RPM cutoff: ", opt_rpm_cutoff))) +
    qc_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          # The first x label is rotated to the left of the panel: widen the margin for long labels
          plot.margin = margin(5.5, 5.5, 5.5, max(5.5, 2.2 * max(nchar(window_qc_df$sample)) - 35), unit = "pt"))

  # B: initiation zones with a used window at each distance, per sample
  qc_coverage_plot <- ggplot(coverage_df %>%
                               mutate(sample = factor(sample, levels = qc_sample_levels),
                                      sample_facet = factor(sample_facet, levels = qc_facet_levels)),
                             aes(x = dist / 1000, y = percent_iz, color = sample, linetype = sample)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.3) +
    scale_color_manual(values = line_colors) +
    scale_linetype_manual(values = line_types) +
    coord_cartesian(ylim = c(0, 100)) +
    facet_wrap(~ sample_facet, nrow = 1, labeller = labeller(.cols = qc_facet_labels), drop = FALSE) +
    xlab("Distance from initiation zone center (kb)") +
    ylab("Initiation zones with a used window (%)") +
    labs(color = caption_lines(paste0("N = ", length(IZ_gr), " initiation zones"),
                               "Line: initiation zones with a used window at each distance"),
         linetype = caption_lines(paste0("N = ", length(IZ_gr), " initiation zones"),
                                  "Line: initiation zones with a used window at each distance")) +
    qc_theme +
    theme(legend.key.width = legend_key_width)

  qc_label_chars <- max(nchar(window_qc_df$sample))
  qc_samples_per_facet <- max(table(window_qc_df$sample_facet[window_qc_df$fate == "used"]))
  qc_width <- 1.2 + length(qc_facet_levels) * max(1.8, 0.6 + 0.35 * qc_samples_per_facet) +
    1.3 + max(0.055 * qc_label_chars, 1.7) + max(0, 0.031 * qc_label_chars - 0.49)
  # The window fates page needs room for the 45-degree sample names
  qc_heights <- c(fates = 2.6 + 0.03 * qc_label_chars, coverage = 2.8)

  # One PDF and one PNG per part
  message("\n[", Sys.time(), "] Saving the window QC plots...")
  qc_plots <- list(fates = qc_fates_plot, coverage = qc_coverage_plot)
  for (part in names(qc_plots)) {
    qc_file <- file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_window_qc_plot.", part))
    ggsave(filename = paste0(qc_file, ".pdf"), plot = qc_plots[[part]],
           width = qc_width, height = qc_heights[[part]], units = "in", limitsize = FALSE)
    ggsave(filename = paste0(qc_file, ".png"), plot = qc_plots[[part]],
           width = qc_width, height = qc_heights[[part]], units = "in", dpi = 600,
           limitsize = FALSE, bg = "white")
  }

} else {

  message("\n[", Sys.time(), "] Skipping the window QC plot: no window lies within the plotted range of the initiation zones.")

}


# ===============================================================================
# Profile difference tests: cluster-based sign-flip permutation test between pairs of samples
# ===============================================================================
# Paired per-bin t statistics; clusters of neighbouring significant bins are tested against
# random sign flips of whole initiation zone curves (Maris & Oostenveld, 2007). Cluster p-values
# already cover the whole distance axis, so they are only adjusted over the pairs (Bonferroni).

message("\n[", Sys.time(), "] Testing the profile differences between samples (cluster-based permutation test)...")

PROFILE_TEST_BIN <- opt_profile_test_bin_size * 1000
PROFILE_PAIRS_PER_FIGURE <- 60
window_widths <- table(partition_df$end - partition_df$start + 1)
if (PROFILE_TEST_BIN < max(as.numeric(names(window_widths)))) {
  warning("\n[", Sys.time(), "] WARNING: --profile_test_bin_size (", opt_profile_test_bin_size,
          " kb) is smaller than the windows of some tables; their empty bins split the clusters.")
}

profile_values_df <- partition_df %>%
  dplyr::filter(RPM >= opt_rpm_cutoff,
                !(duplicate %in% TRUE),
                !is.na(RFD_raw)) %>%
  mutate(bin = floor(dist / PROFILE_TEST_BIN)) %>%
  group_by(sample, sample_facet, break_ID, bin) %>%
  summarise(value = mean(RFD_raw), .groups = "drop") %>%
  mutate(sample = display_name(sample),
         is_rfd = sample %in% rfd_labels)

# Paired t per bin, observed (signs = NULL) or for each row of random signs.
# D0: initiation zone x bin differences (0 if missing); n, ss: per-bin counts and sums of squares.
paired_t <- function(D0, n, ss, signs = NULL) {
  sums <- if (is.null(signs)) matrix(colSums(D0), nrow = 1) else signs %*% D0
  N <- matrix(n, nrow = nrow(sums), ncol = length(n), byrow = TRUE)
  m <- sums / N
  v <- (matrix(ss, nrow = nrow(sums), ncol = length(ss), byrow = TRUE) - N * m^2) / (N - 1)
  t <- m / sqrt(v / N)
  t[!is.finite(t) | N < 3] <- 0
  t
}

# Runs of neighbouring bins beyond the threshold with the same sign
find_clusters <- function(t, thr) {
  s <- ifelse(t > thr, 1L, ifelse(t < -thr, -1L, 0L))
  r <- rle(s)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  keep <- which(r$values != 0)
  tibble(first = starts[keep], last = ends[keep], sign = r$values[keep],
         mass = vapply(keep, function(k) sum(t[starts[k]:ends[k]]), double(1)))
}

# Pairs within the partitions of each facet, and within the RFD tables (once)
profile_sets <- list()
for (facet in unique(profile_values_df$sample_facet)) {
  samples <- intersect(violin_sample_levels,
                       profile_values_df$sample[profile_values_df$sample_facet == facet & !profile_values_df$is_rfd])
  if (length(samples) >= 2) profile_sets[[facet]] <- list(facet = facet, samples = samples)
}
rfd_samples <- intersect(violin_sample_levels, profile_values_df$sample[profile_values_df$is_rfd])
if (length(rfd_samples) >= 2) {
  profile_sets[["RFD"]] <- list(facet = unique(profile_values_df$sample_facet[profile_values_df$is_rfd])[1],
                                samples = rfd_samples)
}

set_label <- function(s) {
  vapply(s, function(x) {
    if (x == "RFD") opt_rfd_exp_type else if (x %in% names(sample_labels)) unname(sample_labels[[x]]) else x
  }, character(1), USE.NAMES = FALSE)
}

set.seed(opt_seed)
profile_clusters_df <- tibble(set = character(), group1 = character(), group2 = character(),
                              n_iz = integer(), cluster = integer(), direction = character(),
                              start_kb = double(), end_kb = double(), n_bins = integer(),
                              mass = double(), p_value = double())
profile_diff_df <- tibble(set = character(), group1 = character(), group2 = character(),
                          bin_start_kb = double(), bin_center_kb = double(), n_iz = integer(),
                          mean_difference = double(), ci_low = double(), ci_high = double(),
                          t = double())

for (set_name in names(profile_sets)) {
  set <- profile_sets[[set_name]]
  set_df <- profile_values_df %>% dplyr::filter(sample_facet == set$facet, sample %in% set$samples)

  for (pair in combn(set$samples, 2, simplify = FALSE)) {

    paired <- dplyr::inner_join(
      set_df %>% dplyr::filter(sample == pair[1]) %>% dplyr::select(break_ID, bin, x = value),
      set_df %>% dplyr::filter(sample == pair[2]) %>% dplyr::select(break_ID, bin, y = value),
      by = c("break_ID", "bin"))
    if (nrow(paired) == 0) next

    bins <- seq(min(paired$bin), max(paired$bin))
    izs <- unique(paired$break_ID)
    D0 <- matrix(0, nrow = length(izs), ncol = length(bins))
    D0[cbind(match(paired$break_ID, izs), match(paired$bin, bins))] <- paired$x - paired$y
    n <- tabulate(match(paired$bin, bins), nbins = length(bins))
    ss <- colSums(D0^2)
    thr <- ifelse(n >= 3, qt(1 - opt_profile_test_threshold / 2, df = pmax(n - 1, 1)), Inf)

    t_obs <- as.vector(paired_t(D0, n, ss))
    clusters <- find_clusters(t_obs, thr)

    # Largest cluster mass of each permutation, with the random signs drawn in chunks to bound memory
    chunks <- diff(c(seq(0, opt_profile_test_permutations - 1, by = 1000), opt_profile_test_permutations))
    null_max <- unlist(lapply(chunks, function(n_perm) {
      signs <- matrix(sample(c(-1, 1), n_perm * length(izs), replace = TRUE), nrow = n_perm)
      apply(paired_t(D0, n, ss, signs), 1, function(t) {
        cl <- find_clusters(t, thr)
        if (nrow(cl) == 0) 0 else max(abs(cl$mass))
      })
    }))

    if (nrow(clusters) > 0) {
      profile_clusters_df <- bind_rows(profile_clusters_df, tibble(
        set = set_name, group1 = pair[1], group2 = pair[2], n_iz = length(izs),
        cluster = seq_len(nrow(clusters)),
        direction = ifelse(clusters$sign > 0, "group1 > group2", "group1 < group2"),
        start_kb = bins[clusters$first] * PROFILE_TEST_BIN / 1000,
        end_kb = (bins[clusters$last] + 1) * PROFILE_TEST_BIN / 1000,
        n_bins = as.integer(clusters$last - clusters$first + 1),
        mass = clusters$mass,
        p_value = (1 + vapply(abs(clusters$mass), function(m) sum(null_max >= m), double(1))) /
                  (opt_profile_test_permutations + 1)))
    } else {
      profile_clusters_df <- bind_rows(profile_clusters_df, tibble(
        set = set_name, group1 = pair[1], group2 = pair[2], n_iz = length(izs),
        cluster = NA_integer_, direction = "no cluster"))
    }

    mean_d <- ifelse(n > 0, colSums(D0) / pmax(n, 1), NA_real_)
    se <- ifelse(n >= 2, sqrt(pmax(ss - n * mean_d^2, 0) / pmax(n - 1, 1) / pmax(n, 1)), NA_real_)
    q <- ifelse(n >= 2, qt(0.975, df = pmax(n - 1, 1)), NA_real_)
    profile_diff_df <- bind_rows(profile_diff_df, tibble(
      set = set_name, group1 = pair[1], group2 = pair[2],
      bin_start_kb = bins * PROFILE_TEST_BIN / 1000,
      bin_center_kb = (bins + 0.5) * PROFILE_TEST_BIN / 1000,
      n_iz = n, mean_difference = mean_d,
      ci_low = mean_d - q * se, ci_high = mean_d + q * se, t = t_obs))
  }
}

profile_clusters_df <- profile_clusters_df %>%
  group_by(set) %>%
  mutate(n_pairs = n_distinct(group1, group2),
         p_adj = pmin(1, p_value * n_pairs)) %>%
  ungroup()

write_tsv(profile_clusters_df,
          file = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_profile_tests.tsv")),
          col_names = TRUE, na = "NA")
write_tsv(profile_diff_df,
          file = file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_profile_differences.tsv")),
          col_names = TRUE, na = "NA")

if (nrow(profile_diff_df) == 0) {

  message("\n[", Sys.time(), "] Skipping the profile difference plot: fewer than two samples to compare in every set.")

} else {

  panel_title <- function(set, g1, g2) paste0(set_label(set), "\n", g1, " -\n", g2)
  profile_diff_df <- profile_diff_df %>%
    mutate(panel = panel_title(set, group1, group2))
  profile_diff_df$panel <- factor(profile_diff_df$panel, levels = unique(profile_diff_df$panel))
  significant_df <- profile_clusters_df %>%
    dplyr::filter(!is.na(p_adj), p_adj <= 0.05) %>%
    mutate(panel = factor(panel_title(set, group1, group2), levels = levels(profile_diff_df$panel)),
           p_label = ifelse(p_adj < 0.001, "p < 0.001", paste0("p = ", signif(p_adj, 2))))

  # One plot per pair, so that every panel has its own axes and axis titles
  plot_pair_difference <- function(panel_name) {
    ggplot(profile_diff_df %>% dplyr::filter(panel == panel_name),
           aes(x = bin_center_kb, y = mean_difference)) +
      geom_rect(data = significant_df %>% dplyr::filter(panel == panel_name),
                aes(xmin = start_kb, xmax = end_kb, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = "#E8C547", alpha = 0.3) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey50", alpha = 0.3) +
      geom_line(linewidth = 0.4) +
      geom_text(data = significant_df %>% dplyr::filter(panel == panel_name),
                aes(x = (start_kb + end_kb) / 2, y = Inf, label = p_label),
                inherit.aes = FALSE, vjust = 1.5, size = 1.8) +
      facet_wrap(~ panel) +
      coord_cartesian(clip = "off") +
      xlab("Distance from initiation zone (kb)") +
      ylab("Mean paired difference (95% CI)") +
      theme_bw(base_family = "Helvetica") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 6),
            axis.text = element_text(size = 6),
            strip.background = element_blank(),
            strip.text = element_text(size = 6),
            aspect.ratio = 1,
            plot.margin = margin(5.5, 10, 5.5, 10, unit = "pt"))
  }

  # One figure per set (partition type, or the RFD tables), in numbered parts of at most
  # PROFILE_PAIRS_PER_FIGURE pairs
  for (set_name in unique(as.character(profile_diff_df$set))) {

    set_panels <- levels(droplevels(profile_diff_df$panel[profile_diff_df$set == set_name]))
    part_panels <- split(set_panels, ceiling(seq_along(set_panels) / PROFILE_PAIRS_PER_FIGURE))
    n_parts <- length(part_panels)
    set_file <- file.path(opt_outdir, paste0(opt_prefix, ".", plot_suffix, "_profile_difference_plot_raw.",
                                             gsub("[^A-Za-z0-9_-]", "_", set_name)))

    for (part in seq_len(n_parts)) {

      panels <- part_panels[[part]]
      n_cols <- min(3, length(panels))
      n_rows <- ceiling(length(panels) / n_cols)

      profile_diff_plot <- ggpubr::ggarrange(plotlist = lapply(panels, plot_pair_difference),
                                             ncol = n_cols, nrow = n_rows) %>%
        ggpubr::annotate_figure(bottom = ggpubr::text_grob(
          caption_lines(paste0("Bins: ", PROFILE_TEST_BIN / 1000, " kb, raw values averaged per initiation zone"),
                        "Line: mean paired difference",
                        "Band: 95% confidence interval",
                        paste0("Shaded: clusters with p <= 0.05 (sign-flip cluster permutation test, ",
                               opt_profile_test_permutations, " permutations, cluster-forming p < ",
                               opt_profile_test_threshold, ")"),
                        "P-values: Bonferroni-adjusted over the pairs"),
          size = 5, family = "Helvetica", hjust = 0, x = 0.01))

      part_file <- if (n_parts > 1) paste0(set_file, ".part", part) else set_file
      message("\n[", Sys.time(), "] Saving profile difference plot (", set_label(set_name),
              if (n_parts > 1) paste0(", part ", part, " of ", n_parts), ")...")
      ggsave(filename = paste0(part_file, ".pdf"),
             plot = profile_diff_plot, width = 0.6 + 2.8 * n_cols, height = 0.7 + 3.1 * n_rows, units = "in",
             limitsize = FALSE)
      ggsave(filename = paste0(part_file, ".png"),
             plot = profile_diff_plot, width = 0.6 + 2.8 * n_cols, height = 0.7 + 3.1 * n_rows, units = "in",
             dpi = 600, limitsize = FALSE, bg = "white")
    }
  }

}


# ===============================================================================
# Scatter plots: RFD (OK-seq) vs partition, one set per OK-seq file
# ===============================================================================

if (HAS_RFD && num_types_with_files > 0) {

  # The partition files of each type are paired by position (the i-th file of each type belongs to
  # the i-th sample)
  label_type <- names(part_files)[lengths(part_files) > 0][1]
  pair_labels <- make.unique(vapply(part_files[[label_type]], sample_name_from_file, character(1), USE.NAMES = FALSE), sep = "_")

  partition_df_flt <- partition_df %>%
    filter(RPM >= opt_rpm_cutoff,
           !is.na(RFD_smooth)) %>%
    mutate(sample = ifelse(is.na(pair_id), sample, pair_labels[pair_id])) %>%
    dplyr::select(interval, sample, sample_facet, RFD_smooth) %>%
    group_by(interval, sample, sample_facet) %>%
    summarise(RFD_smooth_mean = mean(RFD_smooth, na.rm = TRUE), .groups = "drop")

  for (i in seq_along(rfd_labels)) {

    rfd_label <- rfd_labels[i]

    if (!RFD_BIN_SIZES_MATCH[[rfd_label]]) {
      message("\n[", Sys.time(), "] Skipping the scatter plot(s) of ", rfd_label, ": its bin size differs from the partitions', ",
              "so their bins cannot be paired by coordinate.")
      next
    }

    if (!rfd_label %in% partition_df_flt$sample) {
      message("\n[", Sys.time(), "] Skipping the scatter plot(s) of ", rfd_label, ": none of its bins passes the RPM cutoff.")
      next
    }

    # A single OK-seq file keeps the plain file names; several get their sample name added
    scatter_prefix <- if (length(rfd_labels) == 1) opt_prefix else paste0(opt_prefix, ".", rfd_file_ids[i])

    message("\n[", Sys.time(), "] Creating scatter plot(s) (", rfd_label, " vs partitions)...")

    RFD_plot_df <- partition_df_flt %>%
      filter(sample == rfd_label | !sample %in% rfd_labels) %>%
      pivot_wider(names_from = sample, values_from = RFD_smooth_mean) %>%
      pivot_longer(cols = -c(interval, all_of(rfd_label), sample_facet), names_to = "sample", values_to = "Partition") %>%
      dplyr::rename(RFD = all_of(rfd_label)) %>%
      mutate(sample = strip_exp_type(sample, opt_partition_exp_type))

    partition_scatter_plot <- RFD_plot_df %>%
      ggplot(aes(x = RFD, y = Partition)) +
      geom_hex(bins = 100) +
      scale_fill_gradientn(colours = (brewer.pal(n = 9, name = "Blues")[2:8])) +
      geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.5) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.5) +
      geom_smooth(se = FALSE, method = "lm", linewidth = 0.3, color = "red") +
      # add correlation and p-value to the plot
      stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
               method = "spearman", size = 2.5,
               cor.coef.name = c("rho")) +
      xlab(paste(rfd_label, "RFD")) +
      facet_grid(sample ~ sample_facet, scales = "free", labeller = labeller(.cols = sample_labels)) +
      theme_bw(base_family = "Helvetica") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            axis.title = element_text(size = 6),
            axis.text = element_text(size = 6),
            legend.text = element_text(size = 6),
            strip.text = element_text(size = 7),
            strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
            strip.background = element_blank(),
            aspect.ratio = 1)

    scatter_width  <- 1.6 + 1.6 * n_distinct(RFD_plot_df$sample_facet) + 0.06 * max(nchar(RFD_plot_df$sample))
    scatter_height <- 1 + 1.6 * n_distinct(RFD_plot_df$sample)

    message("\n[", Sys.time(), "] Saving scatter plot(s) (", rfd_label, " vs partitions)...")
    ggsave(filename = file.path(opt_outdir, paste0(scatter_prefix, ".scatter_plot.pdf")),
           plot = partition_scatter_plot, width = scatter_width, height = scatter_height, units = "in",
           limitsize = FALSE)

    ggsave(filename = file.path(opt_outdir, paste0(scatter_prefix, ".scatter_plot.png")),
           plot = partition_scatter_plot, width = scatter_width, height = scatter_height, units = "in",
           dpi = 600, limitsize = FALSE)

  }

}

message("\n[", Sys.time(), "] All processing and plotting complete.")
