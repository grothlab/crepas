#!/usr/bin/env Rscript

# ===============================================================================
# repliseq_rtnormalize.R
#
# Written for the grothlab/crepas pipeline by:
#     - Samuel Ruiz-Pérez <samper@cancer.dk>
#     https://github.com/grothlab/crepas/
#
# Description:
#     Calculates and normalizes a replication-timing (E/L ratio) track for one
#     Repli-seq sample from a multi-replicate, multi-phase raw read-count table
#     (chr, start, end, <phase_brep> columns): CPM normalization, optional
#     quantile normalization between replicates, early/late replicate
#     combination (pooled across replicates, paired per biological replicate,
#     or auto-chosen via a correlation heuristic), and loess or rolling-window
#     smoothing applied either to each replicate's coverage or to the ratio.
# ===============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(zoo)
})

# Fewest usable bins a sequence needs before loess is attempted, and the smallest number of
# points its local neighbourhood is allowed to span
MIN_BINS_FOR_LOESS <- 10

# S-phase fractions, earliest first. The order is what gives the replication-timing index its
# meaning; the early/late ratio only ever uses the first and last of these.
PHASE_ORDER <- c("early", "mid", "late")

# ===============================================================================
# Argument parsing
# ===============================================================================

parser <- ArgumentParser()

parser$add_argument("-c", "--counts", action = "store",
                    type = "character",
                    required = TRUE,
                    help = "Raw read-count table (chr, start, end, then one column per replicate), e.g. deepTools multiBamSummary --outRawCounts output; its one header line is always skipped [required]")

parser$add_argument("--phases", action = "store",
                    type = "character",
                    nargs = "+",
                    required = TRUE,
                    help = "Phase ('early' or 'late') of each count column, in the same order as --counts [required]")

parser$add_argument("--breps", action = "store",
                    type = "character",
                    nargs = "+",
                    required = TRUE,
                    help = "Biological replicate identifier of each count column, in the same order as --counts. Used to match early/late replicate pairs for --method paired [required]")

parser$add_argument("-m", "--method", action = "store",
                    default = "pooled",
                    choices = c("auto", "paired", "pooled"),
                    help = "How to combine replicates: 'paired' computes RT per matched early/late biological-replicate pair and averages; 'pooled' averages replicates within each phase first, then computes one RT ratio; 'auto' chooses between the two via a correlation heuristic [default: %(default)s]")

parser$add_argument("-n", "--normalization", action = "store",
                    default = "qn",
                    choices = c("cpm", "qn"),
                    help = "'cpm' only CPM-normalizes each replicate; 'qn' additionally quantile-normalizes between replicates (on top of CPM), within the paired RT columns or within each phase's CPM columns depending on --method [default: %(default)s]")

parser$add_argument("-s", "--smooth", action = "store",
                    default = "loess",
                    choices = c("loess", "roll", "none"),
                    help = "Smoothing strategy [default: %(default)s]")

parser$add_argument("--smooth_stage", action = "store",
                    default = "replicate",
                    choices = c("replicate", "ratio"),
                    help = "Where smoothing enters the calculation: 'replicate' smooths each replicate's normalized coverage track and then takes the ratio (writing out the smoothed per-replicate coverage tracks alongside the ratio); 'ratio' takes the ratio first and smooths that. Only affects the smoothed output; the raw output is always the unsmoothed ratio [default: %(default)s]")

parser$add_argument("--ratio_direction", action = "store",
                    default = "early_over_late",
                    choices = c("early_over_late", "late_over_early"),
                    help = "Orientation of the log2 ratio. 'early_over_late' makes high values early-replicating, which is the usual convention; 'late_over_early' inverts it [default: %(default)s]")

parser$add_argument("--clamp_negative", action = "store",
                    default = TRUE,
                    type = "logical",
                    help = "Whether to floor smoothed replicate coverage at zero before taking the ratio. Coverage cannot be negative, but a loess fit can undershoot below it. Only applies with --smooth_stage replicate [default: %(default)s]")

parser$add_argument("--loess_degree", action = "store",
                    default = 2,
                    type = "integer",
                    choices = c(1, 2),
                    help = "Degree of the local polynomial fitted by loess [default: %(default)s]")

parser$add_argument("--loess_family", action = "store",
                    default = "gaussian",
                    choices = c("gaussian", "symmetric"),
                    help = "Loess fitting family. 'symmetric' reweights to resist outliers, but is undefined on sequences where over half the residuals are exactly zero, in which case this falls back to 'gaussian' [default: %(default)s]")

parser$add_argument("--loess_window", action = "store",
                    default = 500000,
                    type = "double",
                    help = "Width (in bp) of the loess smoothing window. Converted to a per-sequence span, so that the same genomic scale is smoothed over on every sequence regardless of its length [default: %(default)s]")

parser$add_argument("--loess_span", action = "store",
                    type = "double",
                    help = "Loess span given directly, as a fraction of the bins on each sequence. Overrides --loess_window when set. Because it is a fraction rather than a fixed width, a long chromosome is smoothed over a proportionally wider stretch of sequence than a short one [default: unset, i.e. derive the span from --loess_window]")

parser$add_argument("--roll_k", action = "store",
                    default = 5,
                    type = "integer",
                    help = "Rolling-mean window size, in bins [default: %(default)s]")

parser$add_argument("--corr", action = "store",
                    default = "spearman",
                    choices = c("spearman", "pearson"),
                    help = "Correlation method used by the 'auto' replicate-combination heuristic [default: %(default)s]")

parser$add_argument("--pairing_threshold", action = "store",
                    default = 0.02,
                    type = "double",
                    help = "Under --method auto, 'paired' is chosen when (matched-pair avg. correlation - mismatched-pair avg. correlation) is at least this [default: %(default)s]")

parser$add_argument("-g", "--exclude_scaffolds", action = "store",
                    default = TRUE,
                    type = "logical",
                    help = "Whether to exclude scaffolds and alternate sequences from the track. A sequence counts as one when its name begins with 'chrUn', ends with '_random', '_alt' or '_fix', or contains a dot ('.') [default: %(default)s]")

parser$add_argument("--min_bin_reads", action = "store",
                    default = 1,
                    type = "double",
                    help = "Fewest reads a bin needs across the early and late fractions combined to appear in the '.covered' copy of the ratio track. Bins below it carry no timing information: with no reads in either fraction the ratio is log2((0+pseudocount)/(0+pseudocount)) = 0, which is indistinguishable from a genuinely mid-replicating bin. [default: %(default)s]")

parser$add_argument("--pseudocount", action = "store",
                    default = 1,
                    type = "double",
                    help = "Pseudocount added to CPM values before taking log2(early/late) or log10 (for correlation) [default: %(default)s]")

parser$add_argument("--min_points_for_cor", action = "store",
                    default = 50,
                    type = "integer",
                    help = "Minimum number of finite paired observations required to compute a correlation under --method auto; pairs with fewer are treated as NA [default: %(default)s]")

parser$add_argument("-o", "--outdir", action = "store",
                    default = ".",
                    type = "character",
                    help = "Path to output directory [default: %(default)s]")

parser$add_argument("-p", "--prefix", action = "store",
                    default = "",
                    type = "character",
                    help = "Prefix prepended to every output filename, as '<prefix>.<suffix>' [default: none, i.e. just '<suffix>']")

parser$add_argument("--sample_name", action = "store",
                    default = "",
                    type = "character",
                    help = "Name to report this sample under in the MultiQC summary table [default: the value of --prefix]")

opt <- parser$parse_args()

if (!nzchar(opt$sample_name)) {
  opt$sample_name <- opt$prefix
}

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
# Helpers (ported from RepliSeq_process.R, generalized to N replicates)
# ===============================================================================

cpm_vec <- function(x) {
  total <- sum(x, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(x)))
  }
  (x / total) * 1e6
}

safe_cor <- function(x, y, method, min_points) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_points) {
    return(NA_real_)
  }
  cor(x[ok], y[ok], method = method)
}

# Quantile normalization, done here rather than with preprocessCore::normalize.quantiles():
# that function's threaded C implementation aborts with "return code from pthread_create() is
# 22" inside this module's container whatever R_THREADS is set to, so the whole --normalization
# qn path was unusable. Rows holding any non-finite value are left untouched, and ties are
# resolved by interpolating the target distribution at the average rank, which is what
# normalize.quantiles() does too (see the tie handling in preprocessCore's qnorm.c).
qn_matrix <- function(mat) {
  mat <- as.matrix(mat)
  if (ncol(mat) < 2) {
    return(mat)
  }
  ok_rows <- apply(mat, 1, function(x) all(is.finite(x)))
  out <- mat
  if (sum(ok_rows) > 1) {
    sub <- mat[ok_rows, , drop = FALSE]
    # Target distribution: the mean of the sorted columns, i.e. the average k-th order statistic
    target <- rowMeans(apply(sub, 2, sort))
    out[ok_rows, ] <- apply(sub, 2, function(x) {
      approx(seq_along(target), target, xout = rank(x, ties.method = "average"), rule = 2)$y
    })
  }
  out
}

# Replication-timing index: per bin, what proportion of the signal sits in each fraction, then the
# weighted-mean fraction number rescaled to 0..1, where 0 is replicating entirely in the earliest
# fraction, 1 entirely in the latest.
# Bins with no coverage in any fraction have no defined index and are left NA.
rt_index <- function(mat) {
  mat <- as.matrix(mat)
  n_fractions <- ncol(mat)
  out <- rep(NA_real_, nrow(mat))
  total <- rowSums(mat)
  ok <- is.finite(total) & total > 0 & stats::complete.cases(mat)
  if (any(ok)) {
    proportions <- mat[ok, , drop = FALSE] / total[ok]
    weighted_mean_fraction <- as.vector(proportions %*% seq_len(n_fractions))
    out[ok] <- (weighted_mean_fraction - 1) / (n_fractions - 1)
  }
  out
}

# loess(family = "symmetric") derives its robustness weights from the median absolute residual,
# so it fails outright ("NA/NaN/Inf in foreign function call") on any sequence where most bins
# carry no reads in either phase and the RT ratio is therefore exactly 0. Fall back to the
# non-robust fit, and then to no smoothing at all, rather than losing the whole track to one
# degenerate sequence.
smooth_loess_chr <- function(y, x, window, span, degree, family, chr_label) {
  ok <- is.finite(y) & is.finite(x)
  out <- rep(NA_real_, length(y))
  n_ok <- sum(ok)

  if (n_ok < MIN_BINS_FOR_LOESS) {
    message("[", Sys.time(), "] ", chr_label, ": only ", n_ok,
            " usable bins, leaving this sequence unsmoothed.")
    out[ok] <- y[ok]
    return(out)
  }

  # loess takes its span as a fraction of the points. Given a window in bp, convert it to a
  # per-sequence fraction, so that the same genomic scale is smoothed over everywhere; given an
  # explicit span, use it as-is and let the smoothed stretch scale with sequence length. Either
  # way, keep the result between "enough points for loess to fit" and "all of them".
  eff_span <- if (!is.null(span)) {
    span
  } else {
    sequence_span <- diff(range(x[ok]))
    if (sequence_span > 0) window / sequence_span else 1
  }
  eff_span <- min(1, max(eff_span, MIN_BINS_FOR_LOESS / n_ok))

  # Try the requested family first; "gaussian" is always kept as a fallback because "symmetric"
  # is undefined on sequences where over half the residuals are exactly zero.
  for (fam in unique(c(family, "gaussian"))) {
    # A failing fit emits one "all data on boundary of neighborhood" warning per degenerate point
    # on its way to erroring out, which buries the run's real output. Hold each attempt's warnings
    # back and only re-raise them if that attempt is the one whose fit is actually used.
    held_warnings <- character()
    result <- withCallingHandlers(
      tryCatch(
        {
          fit <- loess(y[ok] ~ x[ok], span = eff_span, degree = degree, family = fam,
                       control = loess.control(surface = "direct"))
          predict(fit, x[ok])
        },
        error = function(e) e
      ),
      warning = function(w) {
        held_warnings <<- c(held_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    if (inherits(result, "error") || !any(is.finite(result))) {
      next
    }
    if (fam != family) {
      message("[", Sys.time(), "] ", chr_label,
              ": loess family = \"", family, "\" failed, fell back to \"", fam, "\".")
    }
    for (w in unique(held_warnings)) {
      warning(chr_label, ": loess: ", w, call. = FALSE)
    }
    out[ok] <- result
    return(out)
  }

  message("[", Sys.time(), "] ", chr_label,
          ": loess smoothing failed, leaving this sequence unsmoothed.")
  out[ok] <- y[ok]
  out
}

smooth_roll_chr <- function(y, k) {
  # rollapply returns all-NA once the window is wider than the sequence
  if (length(y) < k) {
    return(y)
  }
  zoo::rollapply(y, width = k, FUN = mean, fill = NA_real_, align = "center", na.rm = TRUE)
}

write_bedgraph <- function(dt, value_col, path) {
  out <- dt[is.finite(get(value_col)), .(chr, start, end, value = get(value_col))]
  setorder(out, chr, start)
  fwrite(out, file = path, sep = "\t", col.names = FALSE)
}

# ===============================================================================
# Read the raw counts table
# ===============================================================================

message("[", Sys.time(), "] Reading counts table: ", opt$counts)
col_ids <- paste(opt$phases, opt$breps, sep = "_")
if (anyDuplicated(col_ids)) {
  stop("Duplicate (phase, brep) combinations in --phases/--breps: ", paste(col_ids[duplicated(col_ids)], collapse = ", "))
}

dt <- fread(opt$counts, header = FALSE, sep = "\t", skip = 1, col.names = c("chr", "start", "end", col_ids))
dt[, start := as.integer(start)]
dt[, end := as.integer(end)]

if (opt$exclude_scaffolds) {
  scaffold_pattern <- "^chrUn|_random$|_alt$|_fix$|\\."
  scaffolds <- unique(dt$chr[grepl(scaffold_pattern, dt$chr)])
  if (length(scaffolds) > 0) {
    message("[", Sys.time(), "] Excluding ", length(scaffolds), " scaffold(s)/alternate sequence(s): ",
            paste(head(scaffolds, 10), collapse = ", "),
            if (length(scaffolds) > 10) ", ..." else "")
    dt <- dt[!grepl(scaffold_pattern, chr)]
  }
  if (nrow(dt) == 0) {
    stop("No bins left after excluding scaffolds. Check the sequence names in the counts table, or pass --exclude_scaffolds FALSE")
  }
}

setorder(dt, chr, start)
dt[, bin_mid := (start + end) / 2]

early_ids <- col_ids[opt$phases == "early"]
late_ids <- col_ids[opt$phases == "late"]
early_breps <- opt$breps[opt$phases == "early"]
late_breps <- opt$breps[opt$phases == "late"]

if (length(early_ids) < 1 || length(late_ids) < 1) {
  stop("At least one 'early' and one 'late' replicate are required")
}

# Fractions actually supplied, earliest first. Intermediate fractions take no part in the
# early/late ratio; they are only used by the replication-timing index below.
present_phases <- PHASE_ORDER[PHASE_ORDER %in% opt$phases]
index_phases <- if (length(present_phases) >= 3) present_phases else character(0)
if (length(index_phases) > 0) {
  message("[", Sys.time(), "] Replication-timing index will be calculated over ",
          length(index_phases), " fractions: ", paste(index_phases, collapse = " < "))
}

# ===============================================================================
# CPM normalization per replicate
# ===============================================================================

for (id in col_ids) {
  dt[, (paste0(id, "_cpm")) := cpm_vec(get(id))]
}

pc <- opt$pseudocount

# ===============================================================================
# Decide how to combine replicates: paired, pooled, or auto (correlation heuristic)
# ===============================================================================

matched_breps <- intersect(early_breps, late_breps)

for (id in col_ids) {
  dt[, (paste0(id, "_logcpm")) := log10(get(paste0(id, "_cpm")) + pc)]
}

correlate_ids <- function(a, b) {
  safe_cor(dt[[paste0(a, "_logcpm")]], dt[[paste0(b, "_logcpm")]], opt$corr, opt$min_points_for_cor)
}

# Every same-phase replicate pair: replicate concordance
concordance_pairs <- list()
for (ids in list(early_ids, late_ids)) {
  if (length(ids) < 2) {
    next
  }
  for (i in seq_len(length(ids) - 1)) {
    for (j in seq(i + 1, length(ids))) {
      concordance_pairs[[length(concordance_pairs) + 1]] <-
        list(a = ids[i], b = ids[j], r = correlate_ids(ids[i], ids[j]))
    }
  }
}
concordance_cors <- vapply(concordance_pairs, function(p) p$r, numeric(1))
min_concordance <- if (length(concordance_cors) > 0) min(concordance_cors, na.rm = TRUE) else NA_real_

# Every early x late pair, split by whether the two share a biological replicate
cross_pairs <- list()
matched_cors <- c()
crossed_cors <- c()
for (e in seq_along(early_ids)) {
  for (l in seq_along(late_ids)) {
    r <- correlate_ids(early_ids[e], late_ids[l])
    is_matched <- early_breps[e] == late_breps[l]
    cross_pairs[[length(cross_pairs) + 1]] <-
      list(a = early_ids[e], b = late_ids[l], r = r, matched = is_matched)
    if (is_matched) {
      matched_cors <- c(matched_cors, r)
    } else {
      crossed_cors <- c(crossed_cors, r)
    }
  }
}
paired_avg <- if (length(matched_cors) > 0) mean(matched_cors, na.rm = TRUE) else NA_real_
crossed_avg <- if (length(crossed_cors) > 0) mean(crossed_cors, na.rm = TRUE) else NA_real_
pairing_delta <- paired_avg - crossed_avg

chosen_method <- opt$method
if (opt$method %in% c("auto", "paired") && length(matched_breps) == 0) {
  if (opt$method == "paired") {
    warning("--method paired requested, but no biological replicate is shared between early and late phases (early breps: ",
            paste(early_breps, collapse = ", "), "; late breps: ", paste(late_breps, collapse = ", "),
            "); falling back to 'pooled'.")
  }
  chosen_method <- "pooled"
} else if (opt$method == "auto") {
  chosen_method <- if (is.finite(pairing_delta) && pairing_delta >= opt$pairing_threshold) "paired" else "pooled"
}

message("[", Sys.time(), "] Replicate-combination method: ", chosen_method)

# ===============================================================================
# Compute the raw RT ratio
# ===============================================================================

# log2 of the early/late ratio, in whichever orientation was asked for
ratio_log2 <- function(early, late) {
  if (opt$ratio_direction == "early_over_late") {
    log2((early + pc) / (late + pc))
  } else {
    log2((late + pc) / (early + pc))
  }
}

# Reduce the per-replicate columns named "<id><suffix>" to one RT value per bin. `apply_qn`
# exists because quantile normalization has to happen once and in the right place: on the
# coverage columns before they are smoothed (--smooth_stage replicate), or on the ratio columns
# afterwards (--smooth_stage ratio).
combine_replicates <- function(suffix, apply_qn) {
  if (chosen_method == "paired") {
    rt_mat <- vapply(matched_breps, function(brep) {
      e_id <- early_ids[early_breps == brep][1]
      l_id <- late_ids[late_breps == brep][1]
      ratio_log2(dt[[paste0(e_id, suffix)]], dt[[paste0(l_id, suffix)]])
    }, numeric(nrow(dt)))
    rt_mat <- as.matrix(rt_mat)
    if (apply_qn && ncol(rt_mat) > 1) {
      rt_mat <- qn_matrix(rt_mat)
    }
    rowMeans(rt_mat, na.rm = TRUE)
  } else {
    early_mat <- as.matrix(dt[, paste0(early_ids, suffix), with = FALSE])
    late_mat <- as.matrix(dt[, paste0(late_ids, suffix), with = FALSE])
    if (apply_qn) {
      early_mat <- qn_matrix(early_mat)
      late_mat <- qn_matrix(late_mat)
    }
    ratio_log2(rowMeans(early_mat, na.rm = TRUE), rowMeans(late_mat, na.rm = TRUE))
  }
}

# Smooth one per-bin vector, per sequence, with whichever smoother was asked for
smooth_by_chr <- function(col_in, col_out) {
  if (opt$smooth == "loess") {
    dt[, (col_out) := smooth_loess_chr(get(col_in), bin_mid, window = opt$loess_window,
                                       span = opt$loess_span, degree = opt$loess_degree,
                                       family = opt$loess_family, chr_label = .BY$chr), by = chr]
  } else if (opt$smooth == "roll") {
    dt[, (col_out) := smooth_roll_chr(get(col_in), k = opt$roll_k), by = chr]
  } else {
    dt[, (col_out) := get(col_in)]
  }
}

qn_wanted <- opt$normalization == "qn"

# One column per fraction, earliest first, each the mean of that fraction's replicates.
fraction_means <- function(suffix, apply_qn) {
  vapply(index_phases, function(phase) {
    phase_ids <- col_ids[opt$phases == phase]
    phase_mat <- as.matrix(dt[, paste0(phase_ids, suffix), with = FALSE])
    if (apply_qn) {
      phase_mat <- qn_matrix(phase_mat)
    }
    rowMeans(phase_mat, na.rm = TRUE)
  }, numeric(nrow(dt)))
}

# The raw track is always the plain unsmoothed ratio, whatever --smooth_stage says
dt[, RT_raw := combine_replicates("_cpm", apply_qn = qn_wanted)]
if (length(index_phases) > 0) {
  dt[, RT_index_raw := rt_index(fraction_means("_cpm", apply_qn = qn_wanted))]
}

# ===============================================================================
# Smoothing
# ===============================================================================

if (opt$smooth == "none") {
  dt[, RT_smooth := RT_raw]
} else {
  smoother_desc <- if (opt$smooth == "loess") {
    paste0("loess ", if (!is.null(opt$loess_span)) paste0("span = ", opt$loess_span)
                     else paste0("window = ", opt$loess_window, " bp"),
           ", degree = ", opt$loess_degree, ", family = ", opt$loess_family)
  } else {
    paste0("rolling mean, k = ", opt$roll_k, " bins")
  }
  message("[", Sys.time(), "] Smoothing per ", opt$smooth_stage, " (", smoother_desc, ")...")

  if (opt$smooth_stage == "ratio") {
    smooth_by_chr("RT_raw", "RT_smooth")
    if (length(index_phases) > 0) {
      smooth_by_chr("RT_index_raw", "RT_index_smooth")
    }
  } else {
    # Smooth each replicate's coverage, then combine. Quantile normalization comes first here,
    # so that replicates are put on a common distribution before being smoothed.
    if (qn_wanted) {
      for (phase in present_phases) {
        phase_ids <- col_ids[opt$phases == phase]
        phase_cpm <- paste0(phase_ids, "_cpm")
        dt[, (paste0(phase_ids, "_norm")) := as.data.table(qn_matrix(dt[, ..phase_cpm]))]
      }
    } else {
      for (id in col_ids) {
        dt[, (paste0(id, "_norm")) := get(paste0(id, "_cpm"))]
      }
    }

    for (id in col_ids) {
      smooth_by_chr(paste0(id, "_norm"), paste0(id, "_sm"))
      if (opt$clamp_negative) {
        # Coverage cannot be negative, but the fit can undershoot below zero
        dt[, (paste0(id, "_sm")) := pmax(get(paste0(id, "_sm")), 0)]
      }
      write_bedgraph(dt, paste0(id, "_sm"), out_path(paste0(id, ".coverage.smooth.bedGraph")))
    }

    # QN has already been applied above, so it must not be applied again to the ratios
    dt[, RT_smooth := combine_replicates("_sm", apply_qn = FALSE)]
    if (length(index_phases) > 0) {
      dt[, RT_index_smooth := rt_index(fraction_means("_sm", apply_qn = FALSE))]
    }
  }
}

if (length(index_phases) > 0 && opt$smooth == "none") {
  dt[, RT_index_smooth := RT_index_raw]
}

# ===============================================================================
# Write outputs
# ===============================================================================

# A bin counts as covered when the fractions that form the ratio have any reads in it at all
el_count_cols <- c(early_ids, late_ids)
dt[, el_reads := rowSums(as.matrix(.SD), na.rm = TRUE), .SDcols = el_count_cols]
covered <- dt[el_reads >= opt$min_bin_reads]
n_dropped <- nrow(dt) - nrow(covered)

write_bedgraph(dt, "RT_raw", out_path("RT.raw.bedGraph"))
write_bedgraph(dt, "RT_smooth", out_path("RT.smooth.bedGraph"))
write_bedgraph(covered, "RT_raw", out_path("RT.raw.covered.bedGraph"))
write_bedgraph(covered, "RT_smooth", out_path("RT.smooth.covered.bedGraph"))
message("[", Sys.time(), "] Covered copy of the ratio track drops ", n_dropped, " of ", nrow(dt),
        " bins with fewer than ", opt$min_bin_reads, " read(s) across early and late")
if (length(index_phases) > 0) {
  write_bedgraph(dt, "RT_index_raw", out_path("RT_index.raw.bedGraph"))
  write_bedgraph(dt, "RT_index_smooth", out_path("RT_index.smooth.bedGraph"))
}

qc_lines <- c(
  paste0("counts:\t", opt$counts),
  paste0("fractions:\t", paste(present_phases, collapse = " < ")),
  paste0("early_replicates:\t", paste(early_ids, collapse = ", ")),
  paste0("late_replicates:\t", paste(late_ids, collapse = ", ")),
  if (length(index_phases) > 0) {
    paste0("rt_index_fractions:\t", paste(index_phases, collapse = " < "))
  } else {
    paste0("rt_index:\tnot calculated (needs at least 3 fractions, found ", length(present_phases), ")")
  },
  paste0("requested_method:\t", opt$method),
  paste0("chosen_method:\t", chosen_method),
  paste0("corr_method:\t", opt$corr),
  vapply(concordance_pairs, function(pair) {
    paste0("cpm_correlation(", pair$a, " vs ", pair$b, "):\t", pair$r)
  }, character(1)),
  paste0("min_replicate_concordance:\t", min_concordance),
  vapply(cross_pairs, function(pair) {
    paste0("cpm_correlation(", pair$a, " vs ", pair$b, "):\t", pair$r)
  }, character(1)),
  paste0("matched_replicate_pairs_avg_corr:\t", paired_avg),
  paste0("mismatched_replicate_pairs_avg_corr:\t", crossed_avg),
  paste0("pairing_delta:\t", pairing_delta),
  paste0("pairing_threshold:\t", opt$pairing_threshold),
  paste0("normalization:\t", opt$normalization),
  paste0("pseudocount:\t", opt$pseudocount),
  paste0("ratio_direction:\t", opt$ratio_direction),
  paste0("smoothing:\t", opt$smooth),
  if (opt$smooth != "none") paste0("smooth_stage:\t", opt$smooth_stage) else "",
  if (opt$smooth != "none" && opt$smooth_stage == "replicate") paste0("clamp_negative:\t", opt$clamp_negative) else "",
  if (opt$smooth == "loess" && is.null(opt$loess_span)) paste0("loess_window:\t", opt$loess_window) else "",
  if (opt$smooth == "loess" && !is.null(opt$loess_span)) paste0("loess_span:\t", opt$loess_span) else "",
  if (opt$smooth == "loess") paste0("loess_degree:\t", opt$loess_degree) else "",
  if (opt$smooth == "loess") paste0("loess_family:\t", opt$loess_family) else "",
  if (opt$smooth == "roll") paste0("roll_k:\t", opt$roll_k) else "",
  paste0("n_bins:\t", nrow(dt)),
  paste0("min_bin_reads:\t", opt$min_bin_reads),
  paste0("n_bins_covered:\t", nrow(covered)),
  paste0("n_bins_uncovered:\t", n_dropped)
)
qc_lines <- qc_lines[nzchar(qc_lines)]
writeLines(qc_lines, con = out_path("qc.txt"))

# One row of the same numbers for MultiQC
format_corr <- function(x) {
  if (is.finite(x)) sprintf("%.4f", x) else ""
}

summary_lines <- c(
  paste(c("Sample", "Replicates", "Min replicate r", "Matched r", "Mismatched r", "Delta", "Bins"), collapse = "\t"),
  paste(c(
    opt$sample_name,
    chosen_method,
    format_corr(min_concordance),
    format_corr(paired_avg),
    format_corr(crossed_avg),
    format_corr(pairing_delta),
    nrow(dt)
  ), collapse = "\t")
)
# Deliberately not prefixed with the sample name: the concatenation downstream goes in filename
# order, and "repliseq_rt_header.txt" has to sort before this file for the header to end up on top.
writeLines(summary_lines, con = file.path(opt$outdir, "rt_summary.tsv"))

message("[", Sys.time(), "] Done!")
