#!/usr/bin/env Rscript
# Author: Nicolas Alcaraz <nicolas.alcaraz@cpr.ku.dk>
# Script for plotting partition plots of SCARseq together with it's stranded Input
# Also plots scatter-correlation plots against OK-seq if provided
# downloaded from https://github.com/grothlab/SCARseq_Pipeline/blob/main/SCAR_partition_plots.R

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

required.libs <- c("GenomicAlignments","GenomicFeatures", "RColorBrewer",
                    "dplyr","ggplot2","ggrepel","ggpubr","ggpmisc","reshape2",
                    "hexbin","latex2exp","argparse")

unavailable.libs <- setdiff(required.libs, rownames(installed.packages()))
if (length(unavailable.libs) > 0) {
    BiocManager::install(unavailable.libs)
}

suppressPackageStartupMessages({
  lapply(required.libs, FUN = function(x) {
    do.call("require", list(x))
  })
})


parser <- ArgumentParser()


parser$add_argument("-a","--scar_partition_file", action = "store",
                    #default = "/maps/projects/dan1/people/ngl887/Projects/KU/H3K9me3_bias/SCAR_seq/MAPPING_mm10/rfd_files/AWEMCM2_SCAR_H3K27ac_442_DMSO_T1_r4_t1_R1.srt.nodup.mm10_PE_smooth_results_w1000_s30_d30_z1.txt.gz",
                    type = "character",
                    help = "Partition file from SCAR-seq [required]")


parser$add_argument("-d","--scarminusinput_partition_file", action = "store",
                    #default = "/maps/projects/dan1/people/ngl887/Projects/KU/H3K9me3_bias/SCAR_seq/MAPPING_mm10/rfd_files/AWEMCM2_SCAR_H3K27ac_442_DMSO_T1_r4_t1_R1.srt.nodup.mm10_PE_SCARminusinput_smooth_results_w1000_s30_d30_z1.txt.gz",
                    type = "character",
                    help = "Partition file from input corrected SCAR-seq [required]")

parser$add_argument("-f","--strandedinput_partition_file", action = "store",
                    #default = "/maps/projects/dan1/people/ngl887/Projects/KU/H3K9me3_bias/SCAR_seq/MAPPING_mm10/rfd_files/AWEMCM2_strandedInput_H3K27ac_442_DMSO_T1_r4_t1_R1.srt.nodup.mm10_PE_smooth_results_w1000_s30_d30_z1.txt.gz",
                    type = "character",
                    help = "Partition file from stranded input  [required]")

parser$add_argument("-k", "--okazaki_file", action = "store",
                    #default = "/maps/projects/dan1/people/ngl887/Projects/KU/rep_chromatin_TC/publication/Wenger_et_al_2023/data/external/Okazaki_mm10_r1_smooth_results_w1000_s30_d30_z1.txt.gz",
                    default = "",
                    type = "character",
                    help = "Okazaki file from OK-seq [optional, required for scatter plot]")

parser$add_argument("-i", "--initiation_zones", action = "store",
                    type = "character",
                    #default = "/maps/projects/dan1/people/ngl887/Projects/KU/rep_chromatin_TC/publication/Wenger_et_al_2023/data/external/Initiation_Zones_mm10_mESC.bed",
                    help = "Bed file with known initiation-zones. Must be provided if no Okazaki partition file is given")

parser$add_argument("-b", "--blacklist", action = "store",
                    #default = "/maps/projects/dan1/people/ngl887/scripts/scar_example/data/mm10.blacklist.bed",
                    type = "character",
                    help = "Blacklist bed file with regions to exclude [optional but strongly recommended]")

parser$add_argument("-n", "--prefix", action = "store",
                    default = "SCAR",
                    type = "character",
                    help = "Prefix for output files and plot title")

parser$add_argument("-o", "--outdir", action = "store",
                    default = "plots",
                    type = "character",
                    help = "Path to output directory for plots")

parser$add_argument("-c", "--cpm_cutoff", action = "store",
                    default = 0.3,
                    type = "double",
                    help = "CPM cutoff for noisy bins [default: 0.3]")

parser$add_argument("-r", "--plot_range", action = "store",
                    default = 100,
                    type = "integer",
                    help = "Distance (KB) surrounding Initation Zones to consider for partition plots [default: 100 KB]")

parser$add_argument("-e", "--exclude_chromosomes", action = "store",
                    default = "chrX,chrY,chrM",
                    type = "character",
                    help = "Chromosomes to exclude from analyses, must be provided comma separated [default: chrX,chrY,chrM]")


opt <- parser$parse_args()


plt.width = 19.3
plt.height = 6.52
part.files <- c()
HAS.SCAR <- FALSE
if (is.null(opt$scar_partition_file)) {
} else if (!file.exists(opt$scar_partition_file)) {
  warning("Partition file not found")
} else {
  part.files <- c(part.files, "SCAR" = opt$scar_partition_file)
  HAS.SCAR <- TRUE
}

HAS.SCARINPUT <- FALSE
if (is.null(opt$scarminusinput_partition_file)) {
} else if (!file.exists(opt$scarminusinput_partition_file)) {
  warning("SCAR input-correct file not found")
} else {
  part.files <- c(part.files, "SCAR_Input_Corrected" = opt$scarminusinput_partition_file)
  HAS.SCARINPUT <- TRUE
}

HAS.INPUT <- FALSE
if (is.null(opt$strandedinput_partition_file)) {
} else if (!file.exists(opt$strandedinput_partition_file)) {
  warning("stranded input partition file not found")
} else {
  part.files <- c(part.files, "strandedInput" = opt$strandedinput_partition_file)
  HAS.INPUT <- TRUE
}

if (length(part.files) == 0) {
  stop("Please provide at least one partition file to create plots")
} else if (length(part.files) == 1) {
  plt.width <- 8.7
} else if (length(part.files) == 2) {
  plt.width <- 13
} else if (length(part.files) == 3) {
  plt.width = 19.3
}

HAS.OKSEQ <- FALSE
if (is.null(opt$okazaki_file)) {
  HAS.OKSEQ <- FALSE
} else if (!file.exists(opt$okazaki_file)) {
  HAS.OKSEQ <- FALSE
} else {
  HAS.OKSEQ <- TRUE
}

HAS.IZ <- FALSE
if (is.null(opt$initiation_zones)) {
  HAS.IZ <- FALSE
} else if (!file.exists(opt$initiation_zones)) {
  HAS.IZ <- FALSE
} else {
  HAS.IZ <- TRUE
}

HAS.BLACKLIST <- FALSE
if (is.null(opt$blacklist)) {
  HAS.BLACKLIST <- FALSE
} else if (!file.exists(opt$blacklist)) {
  HAS.BLACKLIST <- FALSE
} else {
  HAS.BLACKLIST <- TRUE
}


if (!HAS.OKSEQ) {
  if (!HAS.IZ) {
    stop("Please provide either valid Okazaki RFD file or Initiation Zone bed file")
  } else {
    print("Okazaki file not found, using provided Initiation Zones")
    IZ.file <- opt$initiation_zones
  }
} else {
  OK.file <- opt$okazaki_file
}

if (!HAS.BLACKLIST) {
  warning("Blacklist file not provided or non-existant, it's recommended to use a blacklist")
}


PREFIX <- opt$prefix
CPM.cutoff <- opt$cpm_cutoff
KB.RANGE <- opt$plot_range
IZ.LIMITS <- KB.RANGE * 1000
plots.dir <- opt$outdir
if (!dir.exists(plots.dir)) {
  dir.create(plots.dir, recursive = TRUE)
}
chrom.excl <- unique(unlist(strsplit(opt$exclude_chromosomes, ",")))


print(paste("CPM cutoff:",opt$cpm_cutoff))


## =================== Load SCAR partiion file =================================



cls <- c("seqnames","start","end", # Coordinates of bin
         "F","R", # Forward (F) and Reverse (R) raw counts in bin
         "F.cpm","R.cpm", # Forward (F.cpm) and Reverse (R.cpm) CPMs in bin
         "RFD.raw","RFD", # Partition scores computed with raw (RFD.raw) and smoothed (RFD) CPMs
         "RFD.deriv", # Value of the derivative of the partition at this bin
         "score", # not used
         "zero.deriv") # second derivative at this bin



blacklist.gr <- GRanges()
if (HAS.BLACKLIST) {
  blacklist.df <- read.csv(opt$blacklist, header = FALSE, sep = "\t")[,1:3]
  colnames(blacklist.df) <- c("seqnames","start","end")
  blacklist.gr <- makeGRangesFromDataFrame(blacklist.df)
}

if (HAS.OKSEQ) {
  OK.df <- read.csv(OK.file,
                      sep="\t",header=FALSE)
  if (ncol(OK.df) != 12) {
    HAS.OKSEQ <- FALSE
    warning("The okazaki RFD file is not in the correct format, it will not be used")
  }
}

if (!(HAS.IZ) & !(HAS.OKSEQ)) {
  stop("ERROR: To create partition plots, please provide Initiation Zones and/or OK-seq partition file for your species, genome version and cell-type")
}

if (HAS.OKSEQ) {
  colnames(OK.df) <- cls
  OK.df <- dplyr::filter(OK.df, !seqnames %in% chrom.excl)
  OK.df$names <- paste(OK.df$seqnames,OK.df$start,sep=":") %>% paste(.,OK.df$end,sep="-") # genomic identifier
  OK.df$exprs <- OK.df[,"F"] + OK.df[,"R"] # total raw counts
  OK.df$CPM <- OK.df$F.cpm + OK.df$R.cpm # total counts pr. million
  OK.gr <- makeGRangesFromDataFrame(OK.df, keep.extra.columns = TRUE)
  OK.gr <- OK.gr[!overlapsAny(OK.gr, blacklist.gr, minoverlap = 1)]
  OK.gr.tmp <- subset(OK.gr, CPM >= CPM.cutoff)

  ids <- which(OK.gr.tmp$zero.deriv > quantile(OK.gr.tmp$RFD.deriv, probs=0.9,
                                               na.rm=TRUE))
  OK.gr.tmp <- OK.gr.tmp[ids]
  names(ids) <- OK.gr.tmp$names
  OK.BIN.SIZE <- as.integer(trimws(OK.df$end[1])) - as.integer(trimws(OK.df$start[1]))
  OK.red.gr <- GenomicRanges::reduce(OK.gr.tmp, min.gapwidth=OK.BIN.SIZE * 3, with.revmap = TRUE)
  filtered.data <- OK.gr.tmp[sapply(OK.red.gr$revmap, function(x) {
    x[which.max(OK.gr.tmp$RFD.deriv[x])]
  })]

  OK.gr$IZ <- ifelse(OK.gr$names %in% filtered.data$names, TRUE, FALSE)
  OK.gr$sample <- "OK-seq"
}




line.colors <- c("SCAR_Input_Corrected" = "darkgreen",
                 "SCAR" = "darkblue",
                 "strandedInput" = "orange")
line.colors <- line.colors[names(part.files)]
sample.labels <-  c("SCAR_Input_Corrected" =  "SCAR (Input-corrected)",
                    "SCAR" = "SCAR", "strandedInput" = "Stranded input")

sample.labels <- sample.labels[names(part.files)]

if (HAS.OKSEQ) {
 line.colors <- c("OK-seq" = "darkgray", line.colors)
 RFD.cor.melt.df <- c()
 RFD.spear.df <- c()
}


RFD.mean.df <- c()

for (pf in names(part.files)) {
  print(pf)
  SCAR.df <- read.csv(part.files[pf], header = FALSE, sep = "\t")
  if (ncol(SCAR.df) != 12) {
    stop(paste("ERROR:", basename(part.files[pf]), "is not a partition file format"))
  }

  colnames(SCAR.df) <- cls
  SCAR.df <- filter(SCAR.df, !seqnames %in% chrom.excl)
  BIN.SIZE <- as.integer(trimws(SCAR.df$end[1])) - as.integer(trimws(SCAR.df$start[1]))
  print(paste("SCAR bin size", BIN.SIZE))

  if (HAS.OKSEQ) {
    if (OK.BIN.SIZE != BIN.SIZE) {
      HAS.OKSEQ <- FALSE
      warning(paste("Bin size of", basename(part.files[pf]), "doesn't correspond to the same bin size as in the partition file, Okazaki won't be used"))
    }
  }

  SCAR.df$names <- paste(SCAR.df$seqnames,SCAR.df$start,sep=":") %>% paste(.,SCAR.df$end,sep="-")
  SCAR.df$IZ <- NA
  SCAR.df$strand <- "*"
  SCAR.df$exprs <- SCAR.df[,"F"] + SCAR.df[,"R"]
  SCAR.df$CPM <- SCAR.df[,"F.cpm"] + SCAR.df[,"R.cpm"]
  SCAR.df$type <- pf
  SCAR.df$sample <- pf


  if (HAS.OKSEQ) {
    Okazaki.df <- as.data.frame(OK.gr)
    Okazaki.df$type <- pf
    Okazaki.df$sample <- "OK-seq"

    Okazaki.df <- Okazaki.df[,colnames(SCAR.df)]
    RFD.gr <- makeGRangesFromDataFrame(rbind(Okazaki.df, SCAR.df),
                                       keep.extra.columns = TRUE)

    OK.ext.gr <- subset(RFD.gr, IZ & sample == "OK-seq")

  } else {
    RFD.gr <- makeGRangesFromDataFrame(SCAR.df,
                                       keep.extra.columns = TRUE)

    IZ.df <- read.csv(opt$initiation_zones, sep = "\t", header = FALSE)[,1:3]
    colnames(IZ.df) <- c("seqnames","start","end")
    print("loadin IZs")
    IZ.df <-  IZ.df[!IZ.df$seqnames %in% chrom.excl,]
    IZ.gr <- makeGRangesFromDataFrame(IZ.df)

    IZ.gr$sample <- "OK-seq"
    OK.ext.gr <- IZ.gr

  }
  OK.ext.gr$break_start <- start(OK.ext.gr)
  # exclude overlapping initiation zone (within 200000 bp).
  # subset to data within search-space (ok.ext.gr)
  OK.ext.gr <- resize(OK.ext.gr,IZ.LIMITS * 2,fix="center")

  OK.dist <- distanceToNearest(OK.ext.gr)
  OK.ext.gr <- OK.ext.gr[-queryHits(subset(OK.dist,
                                           OK.dist@elementMetadata$distance==0))]
  RFD.gr <- RFD.gr[!overlapsAny(RFD.gr, blacklist.gr, minoverlap = 1)]

  n.izs <- length(OK.ext.gr)
  overlap.pairs <- findOverlaps(OK.ext.gr, RFD.gr)
  RFD.break.all.gr <- RFD.gr[subjectHits(overlap.pairs)]
  RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
  RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]

  tmp.mean.df <- RFD.break.all.gr %>% as.data.frame() %>%
    dplyr::filter(F.cpm >= CPM.cutoff | R.cpm >= CPM.cutoff) %>%
    dplyr::group_by(dist,sample, type) %>%  # rank, enh_active
    # TODO: modified this so column RFD is numeric ("" are replaced with NA too)
    dplyr::mutate(RFD = as.numeric(RFD)) %>%
    dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
      RFD.raw = mean(RFD.raw, na.rm = T),
      RFD = mean(RFD,na.rm = T)) %>% as.data.frame()

  RFD.mean.df <- rbind(RFD.mean.df, tmp.mean.df)

  if (HAS.OKSEQ) {
    print("Computing correlations with OK-seq...")

    RFD.cor.mat <- as.data.frame(RFD.gr) %>%
      filter((F.cpm + R.cpm) >= CPM.cutoff) %>%
      group_by(names, sample) %>%
      # TODO: modified this so column RFD is numeric ("" are replaced with NA too)
      dplyr::mutate(RFD = as.numeric(RFD)) %>%
      summarise(RFD = mean(RFD, na.rm = TRUE)) %>% as.data.frame()

    RFD.cor.df  <- reshape2::dcast(RFD.cor.mat, names ~ sample, value.var = "RFD", fun.aggregate = mean)
    colnames(RFD.cor.df) <- c("names","OKseq", pf)

    cor.comp <- colnames(RFD.cor.df)[-1]
    cor.comp <- cor.comp[order(cor.comp)]
    cor.df <- c()
    for (i in 1:(length(cor.comp)-1)) {
      for (j in (i+1):length(cor.comp)) {
        comp.x <- cor.comp[i]
        comp.y <- cor.comp[j]
        cor.obj <- cor.test(RFD.cor.df[,comp.x], RFD.cor.df[,comp.y], method = "spearman",
                            use = "pairwise.complete.obs")
        tmp.df <- data.frame(cor.x = comp.x, cor.y = comp.y,
                             spearman.cor = cor.obj$estimate,
                             cor.pvalue = cor.obj$p.value,
                             stringsAsFactors = FALSE)
        cor.df <- rbind(cor.df, tmp.df)
      }
    }

    RFD.spear.df <- rbind(RFD.spear.df, cor.df)

    tmp.cor.melt.df <- melt(RFD.cor.df, id.vars = c("names","OKseq"),
                            variable.name = "sample", value.name = "RFD") %>%
      filter(!is.nan(RFD) & !is.nan(OKseq))


    colnames(tmp.cor.melt.df) <- c("names", "RFD", "sample", "Partition")
    tmp.cor.melt.df$sample <- as.character(tmp.cor.melt.df$sample)
    RFD.cor.melt.df <- rbind(RFD.cor.melt.df, tmp.cor.melt.df)
  }
}

if (HAS.OKSEQ) {
  colnames(RFD.spear.df) <- c("OKseq","sample","spearman.cor","cor.pvalue")
  RFD.spear.df$tex.label <- paste0("$\\rho = ", sprintf("%0.4f",RFD.spear.df$spearman.cor), "$")
  RFD.spear.df$tex.label <- paste0("rho == ", sprintf("%0.4f",RFD.spear.df$spearman.cor))
  RFD.cor.melt.df$sample <- factor(RFD.cor.melt.df$sample,
                                   levels =  names(sample.labels))
  RFD.spear.df$sample <- factor(RFD.spear.df$sample,
                                levels = names(sample.labels))
}



RFD.mean.df$sample <- factor(RFD.mean.df$sample,
                             levels = names(line.colors))

RFD.mean.df$type <- factor(RFD.mean.df$type,
                             levels = names(sample.labels))


if (HAS.OKSEQ) {
  smooth.plt <- ggplot(RFD.mean.df,aes(x = dist / 1000, y = RFD, colour=sample)) +
    geom_line(stat = "smooth", linewidth = 1.3, aes(x = dist / 1000, y = RFD,
                                                    colour=sample),
              method = "gam",se=F, inherit.aes = TRUE) +
    xlab(paste0("Distance (kb) from initiation zone center")) +
    ylab("Partition or RFD") +
    ggtitle(PREFIX) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.5) +
    scale_colour_manual(values = line.colors,
                        breaks = "OK-seq") +
    guides(guides(colour = guide_legend(order = 1))) +
    theme_bw(base_size = 20, base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #legend.position = "none",
          legend.title = element_blank(),
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 5 * .pt),
          legend.text = element_text(size = 6 * .pt)) +
    facet_wrap(~type, nrow = 1, labeller = labeller(.cols = sample.labels)) +
    labs(caption = TeX(paste0("$N =",n.izs, "$")),
         family = "Helvetica", size = 5)


  raw.plt <- ggplot(RFD.mean.df,aes(x = dist / 1000, y = RFD, colour = sample)) +
    geom_line(linewidth = 1) +
    xlab(paste0("Distance (kb) from initiation zone center")) +
    ylab("Partition or RFD") +
    ggtitle(PREFIX) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.5) +
    scale_colour_manual(values = line.colors,
                        breaks = "OK-seq") +
    guides(guides(colour = guide_legend(order = 1))) +
    theme_bw(base_size = 20, base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #legend.position = "none",
          legend.title = element_blank(),
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 5 * .pt),
          legend.text = element_text(size = 6 * .pt)) +
    facet_wrap(~type, nrow = 1, labeller = labeller(.cols = sample.labels)) +
    labs(caption = TeX(paste0("$N =",n.izs, "$")),
         family = "Helvetica", size = 5)

  max.part <- max(c(RFD.cor.melt.df$RFD, RFD.cor.melt.df$Partition), na.rm = TRUE)

  partition.scatter.plt <- RFD.cor.melt.df %>%
    ggplot(aes(x = RFD, y = Partition)) +
    geom_hex(bins=120) +
    ggtitle(PREFIX) +
    scale_fill_gradientn(colours = (brewer.pal(n=9,name="Blues")[2:8])) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.5) +
    scale_y_continuous(breaks=c(-0.6, 0, 0.6),
                       limits = c(-max.part,max.part)) +
    scale_x_continuous(breaks=c(-0.6, 0, 0.6),
                       limits =  c(-max.part,max.part)) +
    geom_smooth(se=F,method="lm",size=0.3,colour="red") +
    theme_bw(base_size = 20, base_family = "Helvetica") +
    geom_text(data = RFD.spear.df,
              aes(x = -0.4, y = 0.5,
                  label = tex.label),
              parse = TRUE, size = 5) +
    facet_wrap(~ sample, labeller = labeller(.cols = sample.labels)) +
    theme(panel.spacing = unit(0.3, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 6 * .pt),
          legend.text = element_text(size = 6 * .pt),
          strip.text.x = element_text(size = 5 * .pt),
          strip.text.y = element_text(size = 5 * .pt),
          strip.background = element_rect(colour="black", fill="#e3e3e3"))

  ggsave(filename = file.path(plots.dir,paste0(PREFIX,".scatter_plots.pdf")),
         plot = partition.scatter.plt, width = plt.width, height = plt.height,
         device = cairo_pdf,
         dpi = 300)



} else {
  smooth.plt <- ggplot(RFD.mean.df,aes(x = dist / 1000, y = RFD, colour=sample)) +
    geom_line(stat = "smooth", linewidth = 1.3, aes(x = dist / 1000, y = RFD,
                                                    colour=sample),
              method = "gam",se=F, inherit.aes = TRUE) +
    xlab(paste0("Distance (kb) from initiation zone center")) +
    ylab("Partition or RFD") +
    ggtitle(PREFIX) +
    geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
    scale_colour_manual(values = line.colors) +
    guides(guides(colour = guide_legend(order = 1))) +
    theme_bw(base_size = 20, base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 5 * .pt),
          legend.text = element_text(size = 6 * .pt)) +
    facet_wrap(~type, nrow = 1, labeller = labeller(.cols = sample.labels)) +
    labs(caption = TeX(paste0("$N =",n.izs, "$")),
         family = "Helvetica", size = 5)


  raw.plt <- ggplot(RFD.mean.df,aes(x = dist / 1000, y = RFD, colour = sample)) +
    geom_line(linewidth = 1) +
    xlab(paste0("Distance (kb) from initiation zone center")) +
    ylab("Partition or RFD") +
    ggtitle(PREFIX) +
    geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
    scale_colour_manual(values = line.colors) +
    guides(guides(colour = guide_legend(order = 1))) +
    theme_bw(base_size = 20, base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 5 * .pt),
          legend.text = element_text(size = 6 * .pt)) +
    facet_wrap(~type, nrow = 1, labeller = labeller(.cols = sample.labels)) +
    labs(caption = TeX(paste0("$N =",n.izs, "$")),
         family = "Helvetica", size = 5)


}



ggsave(filename = file.path(plots.dir,paste0(PREFIX,".partition_plots_raw.pdf")),
       plot = raw.plt, width = plt.width, height = plt.height,
       device = cairo_pdf,
       dpi = 300)

ggsave(filename = file.path(plots.dir,paste0(PREFIX,".partition_plots_smoothed.pdf")),
       plot = smooth.plt, width = plt.width, height = plt.height,
       device = cairo_pdf,
       dpi = 300)


