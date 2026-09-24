#!/usr/bin/env Rscript

# ===============================================================================
# partition_or_rfd_plot.R
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

required.libs <- c("tidyverse","GenomicAlignments","GenomicFeatures","argparse")

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
# Argument parsing
# ===============================================================================

parser <- ArgumentParser()

parser$add_argument("-k", "--okseq_rfd_file", action = "store",
                    type = "character",
                    required = TRUE,
                    help = "Okazaki file from OK-seq [required]")

parser$add_argument("-b", "--blacklist", action = "store",
                    type = "character",
                    help = "Blacklist bed file with regions to exclude [optional but strongly recommended]")

parser$add_argument("-s", "--chrom_sizes", action = "store",
                    type = "character",
                    help = "Chromosome sizes file [required]")

parser$add_argument("-n", "--prefix", action = "store",
                    default = "RFD",
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

parser$add_argument("-z", "--zero_deriv_quantile", action = "store",
                    default = 0.9,
                    type = "double",
                    help = "Quantile for zero derivative filtering [default: 0.9]")

parser$add_argument("-e", "--exclude_chromosomes", action = "store",
                    default = "chrX,chrY,chrM",
                    type = "character",
                    help = "Chromosomes to exclude from analyses, must be provided comma separated [default: chrX,chrY,chrM]")

parser$add_argument("-g", "--exclude_scaffolds", action = "store",
                    default = TRUE,
                    type = "logical",
                    help = "Whether to exclude scaffolds from analyses. Chromosomes whose name begins with 'chrUn', ends with '_random', '_alt' or '_fix', or contains a dot ('.') are considered scaffolds [default: TRUE]")

parser$add_argument("-l", "--iz_limits_kb", action = "store",
                    default = 100,
                    type = "integer",
                    help = "Upstream and downstream distance in kb to consider for removing overlapping initiation zones. For example, a value of 100 will remove initiation zones that are within 100 kb upstream or 100 kb downstream of the center of another initiation zone [default: 100]")

opt <- parser$parse_args()

opt_okseq_rfd_file <- opt$okseq_rfd_file
opt_blacklist <- opt$blacklist
opt_chrom_sizes <- opt$chrom_sizes
opt_prefix <- opt$prefix
opt_outdir <- opt$outdir
opt_rpm_cutoff <- opt$rpm_cutoff
opt_zero_deriv_quantile <- opt$zero_deriv_quantile
opt_exclude_chromosomes <- opt$exclude_chromosomes
opt_exclude_scaffolds <- opt$exclude_scaffolds
opt_iz_limits_kb <- opt$iz_limits_kb

# ===============================================================================
# Initialize partition files and flags
# ===============================================================================

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
message("# Extracting initiation zones from OK-seq RFD file...")
message("# ===============================================================================")

ok_base_name <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(opt_okseq_rfd_file))
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
        RPM = fwd_RPM + rev_RPM,               # total counts per million
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

message("\n[", Sys.time(), "] (", ok_base_name, ") Keeping only OK-seq bins with sufficient coverage...")
OK_gr <- subset(OK_gr, RPM >= opt_rpm_cutoff)

message("\n[", Sys.time(), "] (", ok_base_name, ") Keeping only OK-seq bins with a RFD zero derivative above the set quantile threshold (", opt_zero_deriv_quantile, ")...")
OK_gr <- OK_gr[which(OK_gr$zero_deriv > quantile(OK_gr$RFD_deriv,
                                            probs = opt_zero_deriv_quantile,
                                            na.rm = TRUE))]

OK_gr_tmp <- GRanges(seqnames = seqnames(OK_gr),
                        ranges = IRanges(start = start(OK_gr),
                                        end = end(OK_gr)),
                        strand = strand(OK_gr))

message("\n[", Sys.time(), "] (", ok_base_name, ") The number of OK-seq bins after filtering by RPM and RFD zero derivative quantile is ", length(OK_gr), ".")
rtracklayer::export.bed(OK_gr_tmp,
                            con = file.path(opt_outdir, paste0(opt_prefix, ".prefiltered.bed")))
rm(OK_gr_tmp)

# As in Petryk et al. (2018; https://www.science.org/doi/10.1126/science.aau0294#supplementary-materials),
# for initiation zones less than 3 bins apart, the bin with the highest RFD derivative is selected:

message("\n[", Sys.time(), "] (", ok_base_name, ") Merging overlapping/adjacent OK-seq bins with a gap smaller than 3x bin size...")
OK_reduced_gr <- GenomicRanges::reduce(OK_gr,
                                    min.gapwidth = OK_BIN_SIZE * 3,
                                    with.revmap = TRUE)

# For each set of merged bins,
message("\n[", Sys.time(), "] (", ok_base_name, ") Keeping the OK-seq bin with the highest RFD derivative...")
filtered_data <- OK_gr[sapply(OK_reduced_gr$revmap,
                                  function(x) { x[which.max(OK_gr$RFD_deriv[x])] })]

rm(OK_reduced_gr)

# Set these bins as initiation zones
OK_gr$IZ <- ifelse(OK_gr$interval %in% filtered_data$interval, TRUE, FALSE)

message("\n[", Sys.time(), "] (", ok_base_name, ") The number of initiation zones after preprocessing is ", sum(OK_gr$IZ), ".")
IZ_gr <- subset(OK_gr, IZ)

if (length(IZ_gr) == 0) {
  warning("[", Sys.time(), "] (", ok_base_name, ") No initiation zones were found; writing empty initiation-zone BED files.")
}

IZ_gr_tmp <- GRanges(seqnames = seqnames(IZ_gr),
                        ranges = IRanges(start = start(IZ_gr),
                                        end = end(IZ_gr)),
                        strand = strand(IZ_gr))

rtracklayer::export.bed(IZ_gr_tmp,
                            con = file.path(opt_outdir, paste0(opt_prefix, ".init_zones.bed")))
rm(IZ_gr_tmp)

message("\n[", Sys.time(), "] (", ok_base_name, ") Removing overlapping initiation zones (within ",
        opt_iz_limits_kb,
        " kb upstream and ",
        opt_iz_limits_kb,
        " kb downstream of another initiation zone)...")

# Get original start coordinate for each initiation zone
IZ_gr$break_start <- start(IZ_gr)
IZ_gr$break_end <- end(IZ_gr)

# Resizing initiation zones to cover 100 kb upstream and 100 kb downstream
IZ_gr <- resize(IZ_gr, opt_iz_limits_kb * 2 * 1000, fix = "center")

# Resizing can generate bins with negative start positions (out-of-bound), so we trim them
IZ_gr <- trim(IZ_gr)

# Finding the nearest resized initiation zone to each resized initiation zone
IZ_dist <- distanceToNearest(IZ_gr)

# Removing overlapping resized initiation zones
overlapping_hits <- queryHits(subset(IZ_dist, IZ_dist@elementMetadata$distance == 0))
# The following line removes overlapping IZs, and catches the case when there are no overlaps
if (length(overlapping_hits) > 0) {
  IZ_gr <- IZ_gr[-overlapping_hits]
} else {
  message("\n[", Sys.time(), "] (", ok_base_name, ") No overlapping initiation zones found within 100 kb upstream and 100 kb downstream of another initiation zone.")
}

# Prepare IZ for saving
IZ_gr_tmp <- GRanges(seqnames = seqnames(IZ_gr),
                        ranges = IRanges(start = IZ_gr$break_start,
                                        end = IZ_gr$break_end),
                        strand = strand(IZ_gr))

message("\n[", Sys.time(), "] (", ok_base_name, ") The number of initiation zones after removing overlaps within ",
        opt_iz_limits_kb,
        " kb upstream and downstream from the IZ center is: ", length(IZ_gr_tmp), ".")

# Save initiation zones to BED file
rtracklayer::export.bed(IZ_gr_tmp,
                            con = file.path(opt_outdir, paste0(opt_prefix, ".init_zones.rm_overlaps.bed")))

message("\n# ===============================================================================")
message("# Done!")
message("# ===============================================================================")
