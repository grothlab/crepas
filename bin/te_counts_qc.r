#!/usr/bin/env Rscript

# te_counts_qc.r
# Samuel Ruiz-Pérez <samper@cancer.dk>, for grothlab/crepas (https://github.com/grothlab/crepas/)
# PCA code adapted from deseq2_qc.r (nf-core/chipseq: Jose Espinosa-Carrasco, Harshil Patel, Gavin Kelly).
#
# Merges per-sample TEcount/TElocal count tables, normalises them with DESeq2
# and plots a PCA and a sample-to-sample distance heatmap.

library(argparse)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

parser <- ArgumentParser(description="Normalise TEcount/TElocal count tables with DESeq2 and plot a PCA and a sample distance heatmap.")
parser$add_argument("count_files", nargs="+", help="Per-sample TEcount/TElocal count tables (*.cntTable or *.cntTable.gz).")
parser$add_argument("-p", "--outprefix", default="te_counts", help="Output prefix. [default: %(default)s]")
parser$add_argument("-f", "--features", default="all", choices=c("te", "gene", "all"),
                    help="Features used for the transformation, PCA and distances. Size factors always use all features. [default: %(default)s]")
parser$add_argument("-m", "--min_count", type="integer", default=10,
                    help="Drop features with fewer reads than this, summed over all samples. [default: %(default)s]")
parser$add_argument("-n", "--ntop", type="integer", default=500,
                    help="Number of most variable features for the first PCA plot; the second uses all of them. [default: %(default)s]")
parser$add_argument("-r", "--rlog", action="store_true", help="Use rlog instead of the variance stabilising transformation (slow on large tables).")
parser$add_argument("-t", "--fit_type", default="parametric", choices=c("parametric", "local", "mean"),
                    help="DESeq2 dispersion fitType for vst/rlog. [default: %(default)s]")
parser$add_argument("-g", "--group_regex", default="",
                    help="Regex removed from sample names to get the group used to colour the plots, e.g. '_bRep_.*$'. Empty: no groups.")
opt <- parser$parse_args()

count_files <- opt$count_files
if (length(count_files) < 2) {
    stop("At least two count files are needed for a PCA.", call.=FALSE)
}

## Sample names are the file names without the dot-separated suffix shared by all files
samples    <- sub("\\.cntTable(\\.gz)?$", "", basename(count_files))
name_parts <- strsplit(samples, ".", fixed=TRUE)
while (all(lengths(name_parts) > 1) && length(unique(sapply(name_parts, tail, 1))) == 1) {
    name_parts <- lapply(name_parts, head, -1)
}
samples <- make.unique(sapply(name_parts, paste, collapse="."))

read_counts <- function(file) {
    table <- read.delim(file, header=TRUE, quote="\"", check.names=FALSE, colClasses=c("character", "numeric"))
    setNames(table[[2]], table[[1]])
}
counts_list <- lapply(count_files, read_counts)

features <- names(counts_list[[1]])
if (all(sapply(counts_list, function(x) identical(names(x), features)))) {
    count_table <- do.call(cbind, counts_list)
} else {
    features    <- unique(unlist(lapply(counts_list, names)))
    count_table <- sapply(counts_list, function(x) { v <- x[features]; v[is.na(v)] <- 0; v })
}
rownames(count_table) <- features
colnames(count_table) <- samples
count_table <- count_table[, order(samples), drop=FALSE]
samples     <- colnames(count_table)

coldata <- data.frame(sample=samples, row.names=samples)
if (opt$group_regex != "") {
    coldata$group <- sub(opt$group_regex, "", samples)
}

gz <- gzfile(paste0(opt$outprefix, ".counts.tsv.gz"), "w")
write.table(cbind(feature=rownames(count_table), as.data.frame(count_table)), file=gz,
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
close(gz)
saveRDS(DESeqDataSetFromMatrix(countData=count_table, colData=coldata, design=~ 1), file=paste0(opt$outprefix, ".dds.rds"))

count_table <- count_table[rowSums(count_table) >= opt$min_count, , drop=FALSE]

is_te         <- grepl(":", rownames(count_table), fixed=TRUE)
selected      <- switch(opt$features, te=is_te, gene=!is_te, all=rep(TRUE, nrow(count_table)))
feature_label <- switch(opt$features, te="TEs", gene="genes", all="features")
if (sum(selected) < 2) {
    stop(paste0("Fewer than two ", feature_label, " with at least ", opt$min_count, " reads in total."), call.=FALSE)
}

dds <- DESeqDataSetFromMatrix(countData=count_table, colData=coldata, design=~ 1)

## Median-of-ratios size factors need at least one feature with no zero in any sample
dds <- tryCatch(
    estimateSizeFactors(dds),
    error = function(e) {
        message("Median-of-ratios size factors failed (", conditionMessage(e), "); using type='poscounts'.")
        estimateSizeFactors(dds, type="poscounts")
    }
)

normalized_counts <- counts(dds, normalized=TRUE)
gz <- gzfile(paste0(opt$outprefix, ".normalized_counts.tsv.gz"), "w")
write.table(cbind(feature=rownames(normalized_counts), as.data.frame(round(normalized_counts, 4))), file=gz,
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
close(gz)
write.table(data.frame(sample=samples, size_factor=round(sizeFactors(dds), 6)), file=paste0(opt$outprefix, ".size_factors.tsv"),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

## Size factors come from all features; the transformation, PCA and distances only use the selected ones
dds <- dds[selected, ]

## vst() fits the dispersion trend on 1000 features, which keeps TElocal tables fast,
## but needs at least that many; smaller tables fall back to the full transformation
if (!opt$rlog) {
    vst_name <- "vst"
    rld <- tryCatch(
        vst(dds, blind=TRUE, fitType=opt$fit_type),
        error = function(e) {
            message("vst() failed (", conditionMessage(e), "); using varianceStabilizingTransformation().")
            varianceStabilizingTransformation(dds, blind=TRUE, fitType=opt$fit_type)
        }
    )
} else {
    vst_name <- "rlog"
    rld <- rlog(dds, blind=TRUE, fitType=opt$fit_type)
}

saveRDS(rld, file=paste0(opt$outprefix, ".", vst_name, ".rds"))
assay(dds, vst_name) <- assay(rld)

plotPCA_vst <- function (object,  ntop = 500, assay=length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    df         <- cbind( as.data.frame(colData(object)), pca$x)
    ord        <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
    df         <- df[ord,]
    attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
    return(df)
}

pdf(file=paste0(opt$outprefix, ".plots.pdf"), onefile=TRUE, width=7, height=7)

for (n_top_var in c(opt$ntop, Inf)) {
    pca.data      <- plotPCA_vst(dds, assay=vst_name, ntop=n_top_var)
    percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
    plot_subtitle <- ifelse(n_top_var==Inf, paste("All", sum(selected), feature_label), paste("Top", min(n_top_var, sum(selected)), "most variable", feature_label))
    point_aes     <- if ("group" %in% colnames(pca.data)) aes(colour=group) else aes()
    pl <- ggplot(pca.data, aes(PC1, PC2, label=paste0(" ", sample, " "))) +
        geom_point(point_aes, size=3) +
        geom_text(check_overlap=TRUE, vjust=0.5, hjust="inward") +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        labs(title = paste0("First PCs on ", vst_name, "-transformed TE counting data"), subtitle = plot_subtitle) +
        theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
    print(pl)
}

pca.vals           <- pca.data[order(pca.data$sample), c("PC1","PC2")]
colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
pca.vals           <- cbind(sample = rownames(pca.vals), round(pca.vals, 4))
write.table(pca.vals, file = paste0(opt$outprefix, ".pca.vals.txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = TRUE)

sampleDists      <- dist(t(assay(dds, vst_name)))
sampleDistMatrix <- as.matrix(sampleDists)
colors           <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(
    sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors,
    annotation_row=if ("group" %in% colnames(coldata)) coldata[, "group", drop=FALSE] else NA,
    main=paste("Euclidean distance between", vst_name, "of samples,", sum(selected), feature_label)
)
dev.off()

write.table(cbind(sample = rownames(sampleDistMatrix), round(sampleDistMatrix, 4)), file=paste0(opt$outprefix, ".sample.dists.txt"),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

sink(paste0(opt$outprefix, ".R_sessionInfo.log"))
print(sessionInfo())
sink()
