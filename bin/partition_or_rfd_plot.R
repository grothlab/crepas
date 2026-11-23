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

# ### Set user's local R library path to packages
# chooseCRANmirror(ind =  28)
# if (!dir.exists(Sys.getenv("R_LIBS_USER"))) {
#     dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
# }

# ### Install required packages if missing
# if (!require("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager", lib = Sys.getenv("R_LIBS_USER"))
# }

required.libs <- c("tidyverse","GenomicAlignments","GenomicFeatures",
           "RColorBrewer","ggrepel","ggpubr","ggpmisc",
          "hexbin","argparse")

unavailable.libs <- setdiff(required.libs, rownames(installed.packages()))
if (length(unavailable.libs) > 0) {
  message("\n[", Sys.time(), "] Installing missing packages: ", paste(unavailable.libs, collapse = ", "))
  BiocManager::install(unavailable.libs)
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

# Get original start coordinate for each initiation zone
IZ_gr$break_start <- start(IZ_gr)

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


plot_partition <- function(mean_col, sd_col, show_sd) {

  p <- ggplot(partition_mean_df, aes(x = dist / 1000, y = .data[[mean_col]], color = sample, linetype = sample)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = "grey96") +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             fill = "grey96") +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3)

  if (show_sd) {
    p <- p +
      geom_ribbon(aes(ymin = .data[[mean_col]] - .data[[sd_col]],
                      ymax = .data[[mean_col]] + .data[[sd_col]], fill = sample),
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
    labs(caption = paste("N =", length(IZ_gr))) +
    theme_bw(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          legend.key.width = legend_key_width,
          legend.key = element_rect(color = NA, fill = NA),
          strip.background = element_blank(),
          aspect.ratio = 1,
          plot.caption = element_text(size = 6)) +
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

smoothed_plot    <- plot_partition("RFD_smooth_mean", "RFD_smooth_sd", show_sd = FALSE)
smoothed_plot_sd <- plot_partition("RFD_smooth_mean", "RFD_smooth_sd", show_sd = TRUE)

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


# ===============================================================================
# Raw (unsmoothed) partition plot(s)
# ===============================================================================

message("\n[", Sys.time(), "] Creating raw partition plot(s)...")

raw_plot    <- plot_partition("RFD_raw_mean", "RFD_raw_sd", show_sd = FALSE)
raw_plot_sd <- plot_partition("RFD_raw_mean", "RFD_raw_sd", show_sd = TRUE)

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
    labs(caption = paste("N =", length(IZ_gr))) +
    theme_bw(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          legend.key.width = legend_key_width,
          legend.key = element_rect(color = NA, fill = NA),
          strip.background = element_blank(),
          aspect.ratio = 1,
          plot.caption = element_text(size = 6)) +
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
