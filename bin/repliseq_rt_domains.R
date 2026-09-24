#!/usr/bin/env Rscript

# ===============================================================================
# repliseq_rt_domains.R
#
# Written for the grothlab/crepas pipeline by:
#     - Samuel Ruiz-Pérez <samper@cancer.dk>
#     https://github.com/grothlab/crepas/
#
# Description:
#     Segments a Repli-seq log2(early/late) replication-timing track into domains
#     of constant timing with circular binary segmentation (DNAcopy), and calls
#     each domain early or late (optionally with an intermediate class) from its
#     segment mean.
# ===============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(DNAcopy)
  library(ggplot2)
})

parser <- ArgumentParser()

parser$add_argument("-b", "--bedgraph", action = "store",
                    type = "character",
                    required = TRUE,
                    help = "Replication-timing log2(early/late) ratio track, as a 4-column bedGraph [required]")

parser$add_argument("--classification", action = "store",
                    default = "three_way",
                    choices = c("binary", "three_way"),
                    help = "'binary' calls every domain early or late; 'three_way' additionally calls domains whose mean sits within --threshold of zero as intermediate, i.e. replicating throughout S phase [default: %(default)s]")

parser$add_argument("--threshold", action = "store",
                    default = 0.1,
                    type = "double",
                    help = "Half-width of the intermediate band around a segment mean of zero, in log2 ratio units. Under --classification binary only its sign is used, so the value has no effect there [default: %(default)s]")

parser$add_argument("--alpha", action = "store",
                    default = 0.01,
                    type = "double",
                    help = "Significance level for accepting a change point in DNAcopy's circular binary segmentation [default: %(default)s]")

parser$add_argument("--min_width", action = "store",
                    default = 2,
                    type = "integer",
                    help = "Smallest number of bins a segment may contain; DNAcopy accepts 2 to 5 [default: %(default)s]")

parser$add_argument("--undo_sd", action = "store",
                    default = 1,
                    type = "double",
                    help = "Merge adjacent segments whose means differ by fewer than this many standard deviations. Zero or less turns the merging off [default: %(default)s]")

parser$add_argument("--smooth_outliers", action = "store",
                    default = TRUE,
                    type = "logical",
                    help = "Whether to run DNAcopy's smooth.CNA() first, which pulls single-bin outliers towards their neighbours before segmenting [default: %(default)s]")

parser$add_argument("--seed", action = "store",
                    default = 1,
                    type = "integer",
                    help = "Random seed. DNAcopy locates change points with a permutation test, so without a fixed seed the same track segments slightly differently on every run [default: %(default)s]")

parser$add_argument("--sample_name", action = "store",
                    default = "sample",
                    type = "character",
                    help = "Name to report this sample under [default: %(default)s]")

parser$add_argument("--skip_plots", action = "store",
                    default = FALSE,
                    type = "logical",
                    help = "Whether to skip the diagnostic plots [default: %(default)s]")

parser$add_argument("--mqc_max_points", action = "store",
                    default = 5000,
                    type = "integer",
                    help = "Most domains per class to put in the MultiQC box plot. Every value is embedded in the report, so larger classes are thinned to a random sample of this size [default: %(default)s]")

parser$add_argument("--plot_width", action = "store",
                    default = 7,
                    type = "double",
                    help = "Width of each plot, in inches [default: %(default)s]")

parser$add_argument("--plot_height", action = "store",
                    default = 5,
                    type = "double",
                    help = "Height of each plot, in inches [default: %(default)s]")

parser$add_argument("-o", "--outdir", action = "store",
                    default = ".",
                    type = "character",
                    help = "Path to output directory [default: %(default)s]")

parser$add_argument("-p", "--prefix", action = "store",
                    default = "",
                    type = "character",
                    help = "Prefix prepended to every output filename, as '<prefix>.<suffix>' [default: none]")

opt <- parser$parse_args()

if (opt$min_width < 2 || opt$min_width > 5) {
  stop("--min_width must be between 2 and 5, which is what DNAcopy's segment() accepts")
}
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

out_path <- function(suffix) {
  filename <- if (nzchar(opt$prefix)) paste(opt$prefix, suffix, sep = ".") else suffix
  file.path(opt$outdir, filename)
}

# ===============================================================================
# Read the track
# ===============================================================================

message("[", Sys.time(), "] Reading track: ", opt$bedgraph)
bins <- fread(opt$bedgraph, header = FALSE, sep = "\t",
              col.names = c("chr", "start", "end", "value"))
bins <- bins[is.finite(value)]
setorder(bins, chr, start)
if (nrow(bins) == 0) {
  stop("No usable bins in ", opt$bedgraph)
}
message("[", Sys.time(), "] ", nrow(bins), " bins across ", uniqueN(bins$chr), " sequences")

# ===============================================================================
# Circular binary segmentation
# ===============================================================================

cna <- CNA(genomdat = as.numeric(bins$value),
           chrom = as.character(bins$chr),
           maploc = as.integer(bins$start),
           data.type = "logratio",
           sampleid = opt$sample_name)

if (opt$smooth_outliers) {
  message("[", Sys.time(), "] Smoothing single-bin outliers...")
  cna <- smooth.CNA(cna)
}

message("[", Sys.time(), "] Segmenting (alpha = ", opt$alpha, ", min.width = ", opt$min_width,
        ", undo.SD = ", opt$undo_sd, ", seed = ", opt$seed, ")...")
set.seed(opt$seed)
segments <- segment(cna,
                    alpha = opt$alpha,
                    min.width = opt$min_width,
                    undo.splits = if (opt$undo_sd > 0) "sdundo" else "none",
                    undo.SD = opt$undo_sd,
                    verbose = 0)

seg <- as.data.table(segments$output)
setnames(seg, c("chrom", "loc.start", "loc.end", "num.mark", "seg.mean"),
         c("chr", "seg_start", "seg_last_bin_start", "n_bins", "seg_mean"))
seg[, chr := as.character(chr)]

# DNAcopy reports a segment's end as the *start* coordinate of its last bin, because that is the
# maploc it was given. Recover the real end from the bin table, or every domain finishes one bin
# short of where it actually does.
bin_ends <- bins[, .(chr, start, end)]
seg <- merge(seg, bin_ends, by.x = c("chr", "seg_last_bin_start"), by.y = c("chr", "start"),
             all.x = TRUE, sort = FALSE)
setnames(seg, "end", "seg_end")
if (anyNA(seg$seg_end)) {
  stop("Could not map ", sum(is.na(seg$seg_end)), " segment end(s) back to a bin")
}
setorder(seg, chr, seg_start)

# ===============================================================================
# Call each domain
# ===============================================================================

if (opt$classification == "binary") {
  seg[, domain := ifelse(seg_mean >= 0, "early", "late")]
} else {
  seg[, domain := ifelse(seg_mean > opt$threshold, "early",
                  ifelse(seg_mean < -opt$threshold, "late", "mid"))]
}

seg[, length := seg_end - seg_start]
# BED score is a 0..1000 integer, so keep the segment mean itself in a trailing column
seg[, score := pmin(1000L, as.integer(round(abs(seg_mean) * 200)))]

classes <- if (opt$classification == "binary") c("early", "late") else c("early", "mid", "late")

write_bed <- function(dt, path, with_extra) {
  cols <- if (with_extra) {
    dt[, .(chr, seg_start, seg_end, domain, score, strand = ".", seg_mean, n_bins)]
  } else {
    dt[, .(chr, seg_start, seg_end, domain)]
  }
  fwrite(cols, file = path, sep = "\t", col.names = FALSE)
}

write_bed(seg, out_path("RT_domains.bed"), with_extra = TRUE)
for (cls in classes) {
  write_bed(seg[domain == cls], out_path(paste0("RT_domains.", cls, ".bed")), with_extra = FALSE)
}

# ===============================================================================
# QC
# ===============================================================================

# Domain sizes are strongly right-skewed, so compare them with rank-based tests rather than
# anything assuming normality. Kruskal-Wallis asks whether size differs across the classes at all;
# the pairwise Wilcoxon tests say which pairs differ, with Benjamini-Hochberg correction over the
# pairs actually tested.
size_test <- list(method = NA_character_, statistic = NA_real_, p_value = NA_real_, pairs = NULL)
testable <- seg[domain %in% classes]
testable[, domain := factor(domain, levels = classes)]
group_sizes <- testable[, .N, by = domain]
if (nrow(group_sizes[N >= 3]) >= 2) {
  usable <- testable[domain %in% group_sizes[N >= 3, domain]]
  usable[, domain := droplevels(domain)]
  if (nlevels(usable$domain) == 2) {
    tt <- suppressWarnings(wilcox.test(length ~ domain, data = usable))
    size_test$method <- "Wilcoxon rank-sum"
  } else {
    tt <- suppressWarnings(kruskal.test(length ~ domain, data = usable))
    size_test$method <- "Kruskal-Wallis"
  }
  size_test$statistic <- unname(tt$statistic)
  size_test$p_value <- tt$p.value

  combos <- combn(levels(usable$domain), 2, simplify = FALSE)
  if (length(combos) > 1) {
    raw_p <- vapply(combos, function(pair) {
      suppressWarnings(wilcox.test(usable[domain == pair[1], length],
                                   usable[domain == pair[2], length])$p.value)
    }, numeric(1))
    size_test$pairs <- data.table(
      comparison = vapply(combos, function(pair) paste(pair, collapse = " vs "), character(1)),
      p_value = raw_p,
      p_adjusted = p.adjust(raw_p, method = "BH")
    )
  }
}

format_p <- function(p) {
  if (!is.finite(p)) "NA" else if (p < 2.2e-16) "< 2.2e-16" else format.pval(p, digits = 3)
}

total_bp <- sum(seg$length)
qc_lines <- c(
  paste0("track:\t", opt$bedgraph),
  paste0("bins:\t", nrow(bins)),
  paste0("sequences:\t", uniqueN(bins$chr)),
  paste0("classification:\t", opt$classification),
  paste0("threshold:\t", opt$threshold),
  paste0("alpha:\t", opt$alpha),
  paste0("min_width:\t", opt$min_width),
  paste0("undo_sd:\t", opt$undo_sd),
  paste0("seed:\t", opt$seed),
  paste0("smooth_outliers:\t", opt$smooth_outliers),
  paste0("segments:\t", nrow(seg)),
  paste0("median_segment_length:\t", as.integer(median(seg$length))),
  paste0("total_bp:\t", format(total_bp, scientific = FALSE))
)
for (cls in classes) {
  n <- sum(seg$domain == cls)
  bp <- sum(seg$length[seg$domain == cls])
  qc_lines <- c(qc_lines,
                paste0(cls, "_segments:\t", n),
                paste0(cls, "_bp:\t", format(bp, scientific = FALSE)),
                paste0(cls, "_fraction_of_bp:\t", sprintf("%.4f", if (total_bp > 0) bp / total_bp else NA_real_)))
}
qc_lines <- c(qc_lines,
              paste0("size_test:\t", size_test$method),
              paste0("size_test_statistic:\t", if (is.finite(size_test$statistic)) sprintf("%.4f", size_test$statistic) else "NA"),
              paste0("size_test_p_value:\t", format_p(size_test$p_value)))
if (!is.null(size_test$pairs)) {
  for (i in seq_len(nrow(size_test$pairs))) {
    qc_lines <- c(qc_lines,
                  paste0("size_test_pair(", size_test$pairs$comparison[i], ")_p_adjusted:\t",
                         format_p(size_test$pairs$p_adjusted[i])))
  }
}
for (cls in classes) {
  lengths_cls <- seg[domain == cls, length]
  qc_lines <- c(qc_lines,
                paste0(cls, "_median_length:\t",
                       if (length(lengths_cls) > 0) as.integer(median(lengths_cls)) else "NA"))
}
writeLines(qc_lines, con = out_path("RT_domains.qc.txt"))

# ===============================================================================
# Plots
# ===============================================================================

if (!opt$skip_plots) {
  message("[", Sys.time(), "] Drawing plots...")

  seg[, domain := factor(domain, levels = classes)]
  class_colours <- setNames(
    if (length(classes) == 2) c("#2166AC", "#B2182B")
    else c("#2166AC", "#B2ABD2", "#B2182B"),
    classes
  )
  base_theme <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")

  # How much of the genome each class covers, which is the headline number
  bp_by_class <- seg[, .(bp = sum(length), n = .N), by = domain]
  bp_by_class[, fraction := bp / sum(bp)]
  p_bp <- ggplot(bp_by_class, aes(x = domain, y = fraction, fill = domain)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%\n(%d domains)", 100 * fraction, n)), vjust = -0.3, size = 3) +
    scale_fill_manual(values = class_colours) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(title = paste0(opt$sample_name, ": genome covered by each class"),
         subtitle = paste0(nrow(seg), " domains from ", nrow(bins), " bins"),
         x = NULL, y = "Fraction of segmented base pairs") +
    base_theme

  # Domain size by class, with the test result stated on the plot
  size_subtitle <- if (is.finite(size_test$p_value)) {
    paste0(size_test$method, ", p = ", format_p(size_test$p_value))
  } else {
    "not enough domains to test for a size difference"
  }
  p_size <- ggplot(seg, aes(x = domain, y = length, fill = domain)) +
    geom_violin(colour = NA, alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = class_colours) +
    scale_y_log10(labels = function(x) paste0(x / 1000, " kb")) +
    labs(title = paste0(opt$sample_name, ": domain size by class"),
         subtitle = size_subtitle, x = NULL, y = "Domain size") +
    base_theme

  # The values the classes were cut from, as a check that the threshold sits where intended
  p_mean <- ggplot(seg, aes(x = domain, y = seg_mean, fill = domain)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_violin(colour = NA, alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = class_colours) +
    labs(title = paste0(opt$sample_name, ": segment mean by class"),
         subtitle = if (opt$classification == "three_way") {
           paste0("Intermediate band is +/- ", opt$threshold, " around zero")
         } else {
           "Split on the sign of the segment mean"
         },
         x = NULL, y = expression(log[2] * "(early / late)")) +
    base_theme

  pdf(out_path("RT_domains.plots.pdf"), width = opt$plot_width, height = opt$plot_height)
  print(p_bp)
  print(p_size)
  print(p_mean)
  invisible(dev.off())

  # Domain sizes for MultiQC as an interactive box plot. MultiQC's violin plot type takes one
  # value per sample rather than a distribution, so a box plot is what raw sizes can be drawn as.
  set.seed(opt$seed)
  box_data <- lapply(classes, function(cls) {
    values <- seg[domain == cls, length]
    if (length(values) > opt$mqc_max_points) {
      values <- sample(values, opt$mqc_max_points)
    }
    as.integer(values)
  })
  names(box_data) <- classes
  box_data <- box_data[lengths(box_data) > 0]

  if (length(box_data) > 0) {
    json_values <- vapply(names(box_data), function(cls) {
      paste0('"', cls, '": [', paste(box_data[[cls]], collapse = ","), ']')
    }, character(1))
    writeLines(paste0(
      '{"id": "repliseq_domain_size",',
      ' "section_name": "MERGED LIB: Repli-seq domain size by class",',
      ' "description": "size of each replication-timing domain, grouped by the class it was',
      ' called as. ', gsub('"', "'", size_subtitle), '.",',
      ' "plot_type": "box",',
      ' "pconfig": {"id": "repliseq_domain_size_plot",',
      ' "title": "Repli-seq: domain size by class",',
      ' "xlab": "Domain size (bp)"},',
      ' "data": {', paste(json_values, collapse = ","), '}}'
    ), con = out_path("RT_domains_size_mqc.json"))
  }
}

message("[", Sys.time(), "] ", nrow(seg), " segments: ",
        paste(sprintf("%s=%d", classes, vapply(classes, function(c) sum(seg$domain == c), integer(1))),
              collapse = ", "))
message("[", Sys.time(), "] Done!")
