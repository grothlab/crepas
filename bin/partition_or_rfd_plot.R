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
#     Also plots scatter-correlation plots against OK-seq if provided.
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
# TODO: This function handles crepas extension conventions. Perhaps a separate sample_id argument would be better.
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


# ===============================================================================
# Argument parsing
# ===============================================================================

parser <- ArgumentParser()

parser$add_argument("-a","--scar_partition_file",
                    action = "store",
                    type = "character",
                    nargs = '+',
                    help = "Partition file(s) from SCAR-seq. It can be a single file path or a space-separated list of quoted file paths. If multiple files are provided, they will be plotted in the same figure.")

parser$add_argument("-d","--scarminusinput_partition_file", action = "store",
                    type = "character",
                    nargs = '+',
                    help = "Partition file(s) from input corrected SCAR-seq. It can be a single file path or a space-separated list of quoted file paths. If multiple files are provided, they will be plotted in the same figure.")

parser$add_argument("-f","--strandedinput_partition_file", action = "store",
                    type = "character",
                    nargs = '+',
                    help = "Partition file(s) from stranded input. It can be a single file path or a space-separated list of quoted file paths. If multiple files are provided, they will be plotted in the same figure.")

parser$add_argument("-k", "--okseq_rfd_file", action = "store",
                    type = "character",
                    help = "Okazaki file from OK-seq [optional, required for scatter plot]")

parser$add_argument("-i", "--initiation_zones", action = "store",
                    type = "character",
                    required = TRUE,
                    help = "Bed file with known initiation-zones. Must be provided if no Okazaki partition file is given")

parser$add_argument("-b", "--blacklist", action = "store",
                    type = "character",
                    help = "Blacklist bed file with regions to exclude [optional but strongly recommended]")

parser$add_argument("-s", "--chrom_sizes", action = "store",
                    type = "character",
                    help = "Chromosome sizes file [required]")

parser$add_argument("-n", "--prefix", action = "store",
                    default = "SCAR",
                    type = "character",
                    help = "Prefix for output files and plot title [default: SCAR]")

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

opt_scar_partition_file <- opt$scar_partition_file
opt_scarminusinput_partition_file <- opt$scarminusinput_partition_file
opt_strandedinput_partition_file <- opt$strandedinput_partition_file
opt_okseq_rfd_file <- opt$okseq_rfd_file
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

part_files <- list("SCAR" = c(), "SCAR_Input_Corrected" = c(), "strandedInput" = c())
HAS_SCAR <- FALSE
HAS_SCARINPUT <- FALSE
HAS_INPUT <- FALSE

# Check SCAR partition files
if (!is.null(opt_scar_partition_file)) {
  for (file in opt_scar_partition_file) {
    if (!file.exists(file)) {
      stop("Partition file not found: ", file)
    } else {
      part_files[["SCAR"]] <- c(part_files[["SCAR"]], file)
      HAS_SCAR <- TRUE
    }
  }
}

# Check input-corrected SCAR partition files
if (!is.null(opt_scarminusinput_partition_file)) {
  for (file in opt_scarminusinput_partition_file) {
    if (!file.exists(file)) {
      stop("SCAR input-corrected file not found: ", file)
    } else {
      part_files[["SCAR_Input_Corrected"]] <- c(part_files[["SCAR_Input_Corrected"]], file)
      HAS_SCARINPUT <- TRUE
    }
  }
}

# Check stranded input partition files
if (!is.null(opt_strandedinput_partition_file)) {
  for (file in opt_strandedinput_partition_file) {
    if (!file.exists(file)) {
      stop("Stranded input partition file not found: ", file)
    } else {
      part_files[["strandedInput"]] <- c(part_files[["strandedInput"]], file)
      HAS_INPUT <- TRUE
    }
  }
}

# Check OK-seq RFD file
HAS_OKSEQ <- FALSE
if (is.null(opt_okseq_rfd_file)) {
  warning("\n[", Sys.time(), "] OK-seq file not provided.")
} else if (!file.exists(opt_okseq_rfd_file)) {
  warning("\n[", Sys.time(), "] OK-seq file not found: ", opt_okseq_rfd_file)
} else {
  HAS_OKSEQ <- TRUE
}

# Count how many types have at least one file
num_types_with_files <- sum(sapply(part_files, function(x) length(x) > 0))

if (num_types_with_files == 0 && !HAS_OKSEQ) {
  stop("[", Sys.time(), "] ERROR: Please provide at least one partition file or OK-seq file to create plots.")
} else if (num_types_with_files == 0 && HAS_OKSEQ) {
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

# Check initiation zones file
HAS_IZ <- FALSE
if (is.null(opt_initiation_zones)) {
  warning("\n[", Sys.time(), "] Initiation zones file not provided.")
} else if (!file.exists(opt_initiation_zones)) {
  warning("\n[", Sys.time(), "] Initiation zones file not found: ", opt_initiation_zones)
} else {
  HAS_IZ <- TRUE
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
HAS_CHROM_SIZES <- FALSE
if (is.null(opt_chrom_sizes)) {
  stop("[", Sys.time(), "] Chromosome sizes file not provided.")
} else if (!file.exists(opt_chrom_sizes)) {
  stop("\n[", Sys.time(), "] Chromosome sizes file not found: ", opt_chrom_sizes)
} else {
  HAS_CHROM_SIZES <- TRUE
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

if (HAS_CHROM_SIZES) {

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
}


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
  mutate(sample = "OK-seq",
          sample_type = "OK-seq") %>%
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


if (HAS_OKSEQ) {

  message("\n# ===============================================================================")
  message("# STEP 2. Preprocessing of OK-seq (RFD) file...")
  message("# ===============================================================================")


  ok_base_name <- sample_name_from_file(opt_okseq_rfd_file)
  message("\n[", Sys.time(), "] (", ok_base_name, ") Reading OK-seq RFD file...")
  OK_df <- read_tsv(opt_okseq_rfd_file, col_names = cls, show_col_types = FALSE)

  if (ncol(OK_df) != 12) {
    stop("[", Sys.time(), "] (", ok_base_name, ") ERROR: The OK-seq RFD file is not in the correct format (ncol != 12)")
  }

  message("\n[", Sys.time(), "] (", ok_base_name, ") Removing OK-seq bins within excluded chromosomes, outside of chrom_sizes, or with <1 fwd or <1 rev counts...")
  OK_df <- OK_df %>%
    dplyr::filter(seqnames %in% names(chrom_sizes),
                  (fwd_counts >= 1 | rev_counts >= 1)) %>%
    mutate(total_counts = fwd_counts + rev_counts, # total raw counts
           RPM = fwd_RPM + rev_RPM,                 # total counts per million
           sample = "OK-seq",
           sample_type = "OK-seq")

  message("\n[", Sys.time(), "] (", ok_base_name, ") Converting the start positions of the OK-seq bins to 1-based...")
  OK_gr <- makeGRangesFromDataFrame(OK_df,
                                    seqinfo = chrom_sizes,
                                    keep.extra.columns = TRUE,
                                    starts.in.df.are.0based = TRUE)

  # We copy the interval now and not before with dplyr because the start
  # coordinates are now 1-based thanks to starts.in.df.are.0based = TRUE
  OK_gr$interval <- paste0(seqnames(OK_gr), ":", start(OK_gr), "-", end(OK_gr))

  message("\n[", Sys.time(), "] (", ok_base_name, ") Calculating OK-seq bin size...")

  OK_BIN_SIZE <- width(OK_gr)[1]

  message("\n[", Sys.time(), "] (", ok_base_name, ") Removing OK-seq bins that overlap a blacklisted region...")
  if (HAS_BLACKLIST) {
    OK_gr <- OK_gr[!overlapsAny(OK_gr, blacklist_gr, minoverlap = 1)]
  }

  # Finding out which initiation zones overlap which OK-seq bins
  if (opt_only_plot_within_iz) {
    overlap_pairs <- findOverlaps(query = OK_gr, subject = IZ_gr, type = "within")
  } else {
    overlap_pairs <- findOverlaps(query = OK_gr, subject = IZ_gr)
  }

  message("\n[", Sys.time(), "] (", ok_base_name, ") Extracting the OK-seq bins that overlap with an initiation zone...")
  RFD_gr <- OK_gr[queryHits(overlap_pairs)]

  # Adding to each OK-seq bin the interval of the initiation zone with which it overlaps
  RFD_gr$break_ID <- IZ_gr$interval[subjectHits(overlap_pairs)]

  message("\n[", Sys.time(), "] (", ok_base_name, ") Calculating the distance from each OK-seq bin to the initiation zone with which it overlaps...")
  RFD_gr$dist <- start(RFD_gr) - IZ_gr$break_start[subjectHits(overlap_pairs)]

  # Remove temporary variables
  rm(OK_gr, overlap_pairs)

} else {

  message("\n# ===============================================================================")
  message("# STEP 2. SKIPPED (No OK-seq file)")
  message("# ===============================================================================")

}

message("\n# ===============================================================================")
message("# STEP 3. Preprocessing of partition files...")
message("# ===============================================================================")

BIN_SIZES_MATCH <- TRUE

partition_df <- tibble()

for (type in names(part_files)) {

  message("\n[", Sys.time(), "] Processing partition files of type: ", type, "...")

  for (file in part_files[[type]]) {

    base_name <- sample_name_from_file(file)

    message("\n[", Sys.time(), "] (", base_name, ") Reading partition file...")
    SCAR_df <- read_tsv(file, col_names = cls, show_col_types = FALSE)

    if (ncol(SCAR_df) != 12) {
      stop("[", Sys.time(), "] ERROR:", base_name, "is not a partition file format (ncol != 12)")
    }

    message("\n[", Sys.time(), "] (", base_name, ") Removing partition bins within excluded chromosomes, outside of chromosome sizes, or with <1 fwd or <1 rev counts...")
    SCAR_df <- SCAR_df %>%
      dplyr::filter(seqnames %in% names(chrom_sizes),
                  (fwd_counts >= 1 | rev_counts >= 1)) %>%
      mutate(total_counts = fwd_counts + rev_counts, # total raw counts
             RPM = fwd_RPM + rev_RPM,                # total counts pr. million
             sample = base_name,
             sample_type = type,
             sample_facet = type)

    SCAR_gr <- makeGRangesFromDataFrame(SCAR_df,
                                        seqinfo = chrom_sizes,
                                        keep.extra.columns = TRUE,
                                        starts.in.df.are.0based = TRUE)

    # We copy the interval now and not before with dplyr because the start
    # coordinates are now 1-based thanks to starts.in.df.are.0based = TRUE
    SCAR_gr$interval <- paste0(seqnames(SCAR_gr), ":", start(SCAR_gr), "-", end(SCAR_gr))

    message("\n[", Sys.time(), "] (", base_name, ") Removing partition bins that overlap a blacklisted region...")
    SCAR_gr <- SCAR_gr[!overlapsAny(SCAR_gr, blacklist_gr, minoverlap = 1)]

    PART_BIN_SIZE <- width(SCAR_gr)[1]
    message("\n[", Sys.time(), "] (", base_name, ") This partition's bin size is ", PART_BIN_SIZE, " bp.")

    if (HAS_OKSEQ) {
      if (OK_BIN_SIZE != PART_BIN_SIZE) {
        BIN_SIZES_MATCH <- FALSE
        warning("\n[", Sys.time(), "] WARNING: Bin size of ", base_name,
                " (", PART_BIN_SIZE, " bp) is not the same bin size as in the provided OK-seq partition file (", OK_BIN_SIZE, " bp). ",
                "The scatter plot pairs bins by their coordinates and will be skipped; the profile plots are not affected.\n")
      }
    }

  # As in Petryk et al. (2018; https://www-science.org/doi/10.1126/science.aau0294#supplementary-materials):
  message("\n[", Sys.time(), "] (", base_name, ") Calculating partition rates around initiation zones (", opt_plot_range, " kb upstream and ", opt_plot_range, " kb downstream of each IZ) by averaging values within each bin position...")

  # Finding out which initiation zones overlap which partition bins...
  if (opt_only_plot_within_iz) {
    overlap_pairs <- findOverlaps(query = SCAR_gr, subject = IZ_gr, type = "within")
  } else {
    overlap_pairs <- findOverlaps(query = SCAR_gr, subject = IZ_gr)
  }

  message("\n[", Sys.time(), "] (", base_name, ") Extracting the partition bins that overlap with an initiation zone...")
  partition_gr <- SCAR_gr[queryHits(overlap_pairs)]

  # Adding to each partition bin the interval of the initiation zone with which it overlaps
  partition_gr$break_ID <- IZ_gr$interval[subjectHits(overlap_pairs)]

  message("\n[", Sys.time(), "] (", base_name, ") Calculating the distance from each partition bin to the initiation zone with which it overlaps...")
  partition_gr$dist <- start(partition_gr) - IZ_gr$break_start[subjectHits(overlap_pairs)]

  partition_df <- partition_df %>%
    bind_rows(as_tibble(partition_gr))

  message("\n[", Sys.time(), "] (", base_name, ") Finished processing partition file.")

  }

}

if (HAS_OKSEQ) {

  message("\n[", Sys.time(), "] Expanding OK-seq RFD data for each partition (plotting purposes)...")
  # Expand RFD_gr for each unique sample_facet in partition_mean_df
  sample_facets <- unique(partition_df$sample_facet)
  # if sample_facets is empty (there are no SCAR partitions, just OK-seq):
  if (length(sample_facets) == 0 || is.null(sample_facets)) {
    sample_facets <- "OK-seq"
  }

  expanded_RFD <- tidyr::expand_grid(
    as_tibble(RFD_gr),
    sample_facet = sample_facets
  )

  partition_df <- partition_df %>%
    bind_rows(expanded_RFD)

}

message("\n# ===============================================================================")
message("# STEP 4. Calculating mean partition/RFD rates around initiation zones...")
message("# ===============================================================================")
# As in Petryk et al. (2018; https://www-science.org/doi/10.1126/science.aau0294#supplementary-materials):

# Filter by RPM cutoff
partition_mean_df <- partition_df %>%
  dplyr::filter(RPM >= opt_rpm_cutoff)

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
  mutate(sample = gsub("^SCAR-seq_", "", sample))

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
  sample_labels <- c("SCAR_Input_Corrected" = "SCAR (Input-corrected)",
                     "SCAR" = "SCAR", "strandedInput" = "Stranded input")
  sample_labels <- sample_labels[names(part_files)]
} else if (HAS_OKSEQ) {
  # Only OK-seq, no partition files
  sample_labels <- c("OK-seq" = "OK-seq")
}

# set a color for each sample, but set OK-seq to dark gray
sample_names <- unique(partition_mean_df$sample)
sample_colors <- setNames(sample_palette(length(sample_names)), sample_names)

# A legend with many samples takes up more of the figure, so widen it to keep the panel readable
if (length(sample_names) > 15) {
  plot_width <- plot_width + 3
}
# Only set OK-seq to grey if there are partition files
if (num_types_with_files > 0) {
  sample_colors["OK-seq"] <- "grey60"
}
line_colors <- sample_colors


# ===============================================================================
# Smoothed partition plot(s)
# ===============================================================================

message("\n[", Sys.time(), "] Creating smoothed partition plot(s)...")


plot_partition <- function(mean_col, sd_col, show_sd) {

  p <- ggplot(partition_mean_df, aes(x = dist / 1000, y = .data[[mean_col]], color = sample)) +
    geom_rect(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
              fill = "grey95", inherit.aes = FALSE) +
    geom_rect(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
              fill = "grey95", inherit.aes = FALSE) +
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
    scale_fill_manual(values = line_colors) +
    xlab("Distance from initiation zone center (kb)") +
    ylab(ifelse(HAS_OKSEQ, ifelse(num_types_with_files > 0, "Partition or RFD", "RFD"), "Partition")) +
    labs(caption = paste("N =", length(IZ_gr))) +
    theme_bw(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          legend.key = element_rect(color = NA, fill = NA),
          strip.background = element_blank(),
          aspect.ratio = 1,
          plot.caption = element_text(size = 6)) +
    guides(color = guide_legend(order = 1)) +
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

smoothed_gam_plot <- ggplot(partition_mean_df, aes(x = dist / 1000, y = RFD_smooth_mean, color = sample)) +
    geom_rect(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
              fill = "grey95", inherit.aes = FALSE) +
    geom_rect(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
              fill = "grey95", inherit.aes = FALSE) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(stat = "smooth", method = "gam", se = FALSE, linewidth = 0.5) +
    scale_color_manual(values = line_colors) +
    xlab("Distance from initiation zone center (kb)") +
    ylab(ifelse(HAS_OKSEQ, ifelse(num_types_with_files > 0, "Partition or RFD", "RFD"), "Partition")) +
    labs(caption = paste("N =", length(IZ_gr))) +
    theme_bw(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          legend.key = element_rect(color = NA, fill = NA),
          strip.background = element_blank(),
          aspect.ratio = 1,
          plot.caption = element_text(size = 6)) +
    guides(color = guide_legend(order = 1)) +
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
# Scatter plots: RFD (OK-seq) vs partition
# ===============================================================================

if (HAS_OKSEQ && num_types_with_files > 0 && !BIN_SIZES_MATCH) {
  message("\n[", Sys.time(), "] Skipping the scatter plot(s): the OK-seq and partition bin sizes differ, ",
          "so their bins cannot be paired by coordinate.")
}

if (HAS_OKSEQ && num_types_with_files > 0 && BIN_SIZES_MATCH) {

  message("\n[", Sys.time(), "] Creating scatter plot(s) (OK-seq vs partitions)...")

  partition_df_flt <- partition_df %>%
    filter(RPM >= opt_rpm_cutoff,
           !is.na(RFD_smooth)) %>%
    dplyr::select(interval, sample, sample_type, sample_facet, RFD_smooth)

  # To prevent empty facets in the scatter plot...
  # If there are stranded input samples, replace their sample names with corresponding SCAR sample names
  if (HAS_INPUT && HAS_SCAR) {
    # Get the sample names for SCAR and strandedInput in order
    scar_samples <- unique(partition_df_flt$sample[partition_df_flt$sample_type == "SCAR"])
    input_samples <- unique(partition_df_flt$sample[partition_df_flt$sample_type == "strandedInput"])

    # Create a mapping from input sample names to SCAR sample names
    if (length(scar_samples) == length(input_samples)) {
      sample_mapping <- setNames(scar_samples, input_samples)

      # Replace strandedInput sample names
      partition_df_flt <- partition_df_flt %>%
        mutate(sample = ifelse(sample_type == "strandedInput",
                               sample_mapping[sample],
                               sample))
    }
  }

  RFD_plot_df <- partition_df_flt %>%
    group_by(interval, sample, sample_facet) %>%
    summarise(RFD_smooth_mean = mean(RFD_smooth, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = RFD_smooth_mean) %>%
    pivot_longer(cols = -c(interval, "OK-seq", sample_facet), names_to = "sample", values_to = "Partition") %>%
    dplyr::rename(RFD = "OK-seq") %>%
    mutate(sample = gsub("^SCAR-seq_", "", sample))

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
    theme_bw(base_size = 20, base_family = "Helvetica") +
    facet_grid(sample ~ sample_facet, scales = "free", labeller = labeller(.cols = sample_labels)) +
    theme_bw(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          strip.background = element_blank(),
          aspect.ratio = 1)

  message("\n[", Sys.time(), "] Saving scatter plot(s) (OK-seq vs partitions)...")
  ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".scatter_plot.pdf")),
         plot = partition_scatter_plot, width = plot_width, height = 6.52, units = "in")

  ggsave(filename = file.path(opt_outdir, paste0(opt_prefix, ".scatter_plot.png")),
         plot = partition_scatter_plot, width = plot_width, height = 6.52, units = "in",
         dpi = 600)

}

message("\n[", Sys.time(), "] All processing and plotting complete.")
