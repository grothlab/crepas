#!/usr/bin/env Rscript

# ===============================================================================
# repliseq_classify_genes.R
#
# Written for the grothlab/crepas pipeline by:
#     - Samuel Ruiz-Pérez <samper@cancer.dk>
#     https://github.com/grothlab/crepas/
#
# Description:
#     Classifies each gene by the S-phase fraction with the highest read density
#     over its gene body, and reports how far that agrees with the call made from
#     the early/late ratio alone.
# ===============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(ggplot2)
})

PHASE_ORDER <- c("early", "mid", "late")

parser <- ArgumentParser()

parser$add_argument("-c", "--counts", action = "store",
                    type = "character",
                    required = TRUE,
                    help = "featureCounts output over gene bodies: Geneid, Chr, Start, End, Strand, Length, then one column per replicate [required]")

parser$add_argument("--phases", action = "store",
                    type = "character",
                    nargs = "+",
                    required = TRUE,
                    help = "Fraction ('early', 'mid' or 'late') of each count column, in the same order as the count columns of --counts [required]")

parser$add_argument("--breps", action = "store",
                    type = "character",
                    nargs = "+",
                    required = TRUE,
                    help = "Biological replicate identifier of each count column, in the same order as the count columns of --counts [required]")

parser$add_argument("--min_gene_reads", action = "store",
                    default = 20,
                    type = "double",
                    help = "Fewest reads a gene needs summed over every fraction before it is classified. Genes below it are reported as 'unclassified': with only a handful of reads the fraction with the highest density is decided by sampling noise [default: %(default)s]")

parser$add_argument("--min_margin", action = "store",
                    default = 0,
                    type = "double",
                    help = "Fewest CPM the winning fraction must exceed the runner-up by. Zero accepts any winner [default: %(default)s]")

parser$add_argument("--pseudocount", action = "store",
                    default = 1,
                    type = "double",
                    help = "Pseudocount added before taking the log2 early/late ratio used for the concordance check [default: %(default)s]")

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
                    help = "Most values per class to put in the MultiQC box plot. Every value is embedded in the report, so larger classes are thinned to a random sample of this size [default: %(default)s]")

parser$add_argument("--seed", action = "store",
                    default = 1,
                    type = "integer",
                    help = "Random seed used when thinning values for the MultiQC box plot [default: %(default)s]")

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

if (length(opt$phases) != length(opt$breps)) {
  stop("--phases and --breps must have the same length")
}
if (!all(opt$phases %in% PHASE_ORDER)) {
  stop("--phases values must be one of: ", paste(PHASE_ORDER, collapse = ", "))
}
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

out_path <- function(suffix) {
  filename <- if (nzchar(opt$prefix)) paste(opt$prefix, suffix, sep = ".") else suffix
  file.path(opt$outdir, filename)
}

# ===============================================================================
# Read the counts
# ===============================================================================

message("[", Sys.time(), "] Reading gene counts: ", opt$counts)
# featureCounts writes a "# Program:..." provenance line above the header
counts <- fread(opt$counts, header = TRUE, sep = "\t", skip = "Geneid")
annotation_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
if (!all(annotation_cols %in% names(counts))) {
  stop("Expected featureCounts columns not found: ", paste(setdiff(annotation_cols, names(counts)), collapse = ", "))
}

sample_cols <- setdiff(names(counts), annotation_cols)
if (length(sample_cols) != length(opt$phases)) {
  stop("--phases has ", length(opt$phases), " entries but --counts has ", length(sample_cols),
       " count columns (", paste(sample_cols, collapse = ", "), ")")
}

col_ids <- paste(opt$phases, opt$breps, sep = "_")
setnames(counts, sample_cols, col_ids)
present_phases <- PHASE_ORDER[PHASE_ORDER %in% opt$phases]
message("[", Sys.time(), "] ", nrow(counts), " genes, fractions: ",
        paste(present_phases, collapse = " < "))

# A gene body can be reported as several intervals when featureCounts merges overlapping features,
# so keep the outermost coordinates and the first sequence name for the BED outputs.
first_of <- function(x) sub(";.*$", "", as.character(x))
counts[, `:=`(chr = first_of(Chr),
              gene_start = vapply(strsplit(as.character(Start), ";"), function(v) min(as.numeric(v)), numeric(1)),
              gene_end = vapply(strsplit(as.character(End), ";"), function(v) max(as.numeric(v)), numeric(1)),
              strand = first_of(Strand))]

# ===============================================================================
# Normalize and average within each fraction
# ===============================================================================

# Each column is scaled by its own library size, so fractions sequenced to different depths are
# comparable. Gene length is the same for every fraction of a gene and so cancels out of the
# comparison, but it is reported as a density anyway to keep the numbers interpretable.
for (id in col_ids) {
  total <- sum(counts[[id]], na.rm = TRUE)
  counts[, (paste0(id, "_cpm")) := if (total > 0) (get(id) / total) * 1e6 else NA_real_]
}

for (phase in present_phases) {
  phase_ids <- col_ids[opt$phases == phase]
  phase_mat <- as.matrix(counts[, paste0(phase_ids, "_cpm"), with = FALSE])
  counts[, (paste0(phase, "_cpm")) := rowMeans(phase_mat, na.rm = TRUE)]
  counts[, (paste0(phase, "_density")) := get(paste0(phase, "_cpm")) / (Length / 1000)]
}

phase_cpm_cols <- paste0(present_phases, "_cpm")
phase_mat <- as.matrix(counts[, ..phase_cpm_cols])

counts[, total_reads := rowSums(as.matrix(.SD), na.rm = TRUE), .SDcols = col_ids]

# ===============================================================================
# Classify by the fraction with the highest density
# ===============================================================================

winner_idx <- max.col(phase_mat, ties.method = "first")
sorted_desc <- t(apply(phase_mat, 1, sort, decreasing = TRUE))
margin <- if (ncol(phase_mat) >= 2) sorted_desc[, 1] - sorted_desc[, 2] else rep(Inf, nrow(phase_mat))

counts[, rt_class := present_phases[winner_idx]]
counts[, margin_cpm := margin]
counts[total_reads < opt$min_gene_reads | margin_cpm < opt$min_margin | !is.finite(margin_cpm),
       rt_class := "unclassified"]

# The traditional call, from the early and late fractions only, for the concordance check
counts[, el_log2_ratio := log2((early_cpm + opt$pseudocount) / (late_cpm + opt$pseudocount))]
counts[, el_class := ifelse(el_log2_ratio >= 0, "early", "late")]
counts[rt_class == "unclassified", el_class := NA_character_]

# Reported is the median log2(early/late) ratio of each class, which is expected to be
# positive for early, negative for late, and near zero for the intermediate class.
class_el_stats <- function(cls) {
  values <- counts[rt_class == cls, el_log2_ratio]
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(list(n = 0L, median = NA_real_, median_abs = NA_real_))
  }
  list(n = length(values), median = median(values), median_abs = median(abs(values)))
}

# ===============================================================================
# Write outputs
# ===============================================================================

out_cols <- c("Geneid", "chr", "gene_start", "gene_end", "strand", "Length", "total_reads",
              paste0(present_phases, "_cpm"), paste0(present_phases, "_density"),
              "rt_class", "margin_cpm", "el_log2_ratio", "el_class")
fwrite(counts[, ..out_cols], file = out_path("gene_RT_class.tsv"), sep = "\t")

classes <- c(present_phases, "unclassified")
for (cls in classes) {
  sub <- counts[rt_class == cls]
  fwrite(sub[, .(chr, start = as.integer(gene_start - 1), end = as.integer(gene_end),
                 name = Geneid, score = 0L, strand)],
         file = out_path(paste0("gene_RT_class.", cls, ".bed")), sep = "\t", col.names = FALSE)
}

qc_lines <- c(
  paste0("counts:\t", opt$counts),
  paste0("sample:\t", opt$sample_name),
  paste0("fractions:\t", paste(present_phases, collapse = " < ")),
  paste0("genes:\t", nrow(counts)),
  paste0("min_gene_reads:\t", opt$min_gene_reads),
  paste0("min_margin:\t", opt$min_margin)
)
for (cls in classes) {
  n <- sum(counts$rt_class == cls)
  qc_lines <- c(qc_lines,
                paste0(cls, "_genes:\t", n),
                paste0(cls, "_fraction_of_genes:\t", sprintf("%.4f", n / nrow(counts))))
}
for (cls in present_phases) {
  st <- class_el_stats(cls)
  qc_lines <- c(qc_lines,
                paste0(cls, "_median_el_log2_ratio:\t", if (is.finite(st$median)) sprintf("%.4f", st$median) else "NA"),
                paste0(cls, "_median_abs_el_log2_ratio:\t", if (is.finite(st$median_abs)) sprintf("%.4f", st$median_abs) else "NA"))
}
mid_stats <- class_el_stats("mid")
end_abs <- c(class_el_stats("early")$median_abs, class_el_stats("late")$median_abs)
mid_is_intermediate <- is.finite(mid_stats$median_abs) && all(is.finite(end_abs)) &&
  mid_stats$median_abs < min(end_abs)
qc_lines <- c(qc_lines,
              paste0("mid_closer_to_el_boundary_than_early_and_late:\t",
                     if ("mid" %in% present_phases) mid_is_intermediate else "NA (no intermediate fraction)"))
writeLines(qc_lines, con = out_path("gene_RT_class.qc.txt"))

# Gene counts per class for MultiQC, which stacks the numeric columns into one bar per sample.
summary_lines <- c(
  paste(c("Sample", classes), collapse = "\t"),
  paste(c(opt$sample_name,
          vapply(classes, function(c) sum(counts$rt_class == c), integer(1))), collapse = "\t")
)
writeLines(summary_lines, con = file.path(opt$outdir, "rt_gene_class_counts.tsv"))

message("[", Sys.time(), "] ",
        paste(sprintf("%s=%d", classes, vapply(classes, function(c) sum(counts$rt_class == c), integer(1))),
              collapse = ", "))
for (cls in present_phases) {
  st <- class_el_stats(cls)
  message("[", Sys.time(), "] ", cls, ": median log2(early/late) = ",
          if (is.finite(st$median)) sprintf("%+.3f", st$median) else "NA",
          " over ", st$n, " genes")
}
# ===============================================================================
# Plots
# ===============================================================================

if (!opt$skip_plots) {
  message("[", Sys.time(), "] Drawing plots...")

  plot_classes <- c(present_phases, "unclassified")
  # Colour the fractions along the S-phase axis and leave the unclassified genes grey
  class_colours <- setNames(
    c(colorRampPalette(c("#2166AC", "#B2ABD2", "#B2182B"))(length(present_phases)), "grey70"),
    plot_classes
  )
  counts[, rt_class := factor(rt_class, levels = plot_classes)]

  base_theme <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")

  # How many genes fall in each class
  class_counts <- counts[, .N, by = rt_class]
  p_counts <- ggplot(class_counts, aes(x = rt_class, y = N, fill = rt_class)) +
    geom_col() +
    geom_text(aes(label = N), vjust = -0.3, size = 3) +
    scale_fill_manual(values = class_colours) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = paste0(opt$sample_name, ": genes per replication-timing class"),
         subtitle = "Assigned to the S-phase fraction with the highest read density over the gene body",
         x = NULL, y = "Genes") +
    base_theme

  # Where each class sits relative to the early/late boundary. This is the plot behind the
  # numbers in the QC file: early above zero, late below, the intermediate class close to it.
  p_el <- ggplot(counts[rt_class != "unclassified" & is.finite(el_log2_ratio)],
                 aes(x = rt_class, y = el_log2_ratio, fill = rt_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_violin(colour = NA, alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = class_colours) +
    labs(title = paste0(opt$sample_name, ": early/late ratio by class"),
         #subtitle = "The two calls cannot disagree on early and late, so what matters is that the intermediate class sits near zero",
         x = NULL, y = expression(log[2] * "(early / late)")) +
    base_theme

  # The same classification seen in the coverage the call was made from
  p_scatter <- ggplot(counts[is.finite(early_cpm) & is.finite(late_cpm)],
                      aes(x = early_cpm + 1, y = late_cpm + 1, colour = rt_class)) +
    geom_point(size = 0.4, alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_manual(values = class_colours) +
    labs(title = paste0(opt$sample_name, ": early against late coverage"),
         subtitle = "Genes on the dashed line have equal early and late signal",
         x = "Early CPM + 1", y = "Late CPM + 1", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

  pdf(out_path("gene_RT_class.plots.pdf"), width = opt$plot_width, height = opt$plot_height)
  print(p_counts)
  print(p_el)
  print(p_scatter)
  invisible(dev.off())

  # The same distribution for MultiQC as an interactive box plot.
  set.seed(opt$seed)
  box_data <- lapply(present_phases, function(cls) {
    values <- counts[rt_class == cls & is.finite(el_log2_ratio), el_log2_ratio]
    if (length(values) > opt$mqc_max_points) {
      values <- sample(values, opt$mqc_max_points)
    }
    round(values, 4)
  })
  names(box_data) <- present_phases
  box_data <- box_data[lengths(box_data) > 0]

  if (length(box_data) > 0) {
    json_values <- vapply(names(box_data), function(cls) {
      paste0('"', cls, '": [', paste(box_data[[cls]], collapse = ","), ']')
    }, character(1))
    writeLines(paste0(
      '{"id": "repliseq_gene_class_el",',
      ' "section_name": "MERGED LIB: Repli-seq early/late ratio by gene class",',
      ' "description": "distribution of the log2(early/late) ratio of the genes assigned to each',
      ' S-phase fraction. Early is expected above zero, late below, and the intermediate class',
      ' close to it. Agreement between the two calls is not shown because it cannot be measured:',
      ' a gene assigned to the early fraction has, by construction, more signal in early than in',
      ' late, which is the only thing the early/late call evaluates.",',
      ' "plot_type": "box",',
      ' "pconfig": {"id": "repliseq_gene_class_el_plot",',
      ' "title": "Repli-seq: early/late ratio by gene class",',
      ' "xlab": "log2(early / late)"},',
      ' "data": {', paste(json_values, collapse = ","), '}}'
    ), con = out_path("gene_RT_class_el_mqc.json"))
  }
}

message("[", Sys.time(), "] Done!")
