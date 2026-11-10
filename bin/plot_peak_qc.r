#!/usr/bin/env Rscript

# ===============================================================================
# plot_peak_qc.r
#
# Originally created by:
#     - Jose Espinosa-Carrasco <https://github.com/JoseEspinosa>
#     - Harshil Patel <https://github.com/drpatelh>
# Source:
#     https://github.com/nf-core/chipseq/blob/76e2382b6d443db4dc2396e6831d1243256d80b0/bin/plot_macs3_qc.r
#
# Adapted for the grothlab/crepas pipeline by:
#     - Samuel Ruiz-Pérez <samper@cancer.dk>
#     https://github.com/grothlab/crepas/
#
# Description:
#     Plots the peak count, peak length and the score columns given with --columns
#     of the peak files of any peak caller, and writes their summary statistics.
# ===============================================================================

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(scales)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--peak_files"), type="character", default=NULL, help="Comma-separated list of peak files.", metavar="path"),
                    make_option(c("-s", "--sample_ids"), type="character", default=NULL, help="Comma-separated list of sample ids associated with peak files. Must be unique and in same order as peaks files input.", metavar="string"),
                    make_option(c("-c", "--columns"), type="character", default='', help="Comma-separated score columns to plot, each as column:name[:scale], where column is the 1-based column number and scale is log2 or log10 (log axis) or neglog10 (plot -log10 of the values). Peak count and length are always plotted.", metavar="string"),
                    make_option(c("-o", "--outdir"), type="character", default='./', help="Output directory", metavar="path"),
                    make_option(c("-p", "--outprefix"), type="character", default='peak_qc', help="Output prefix", metavar="string"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$peak_files)){
    print_help(opt_parser)
    stop("At least one peak file must be supplied", call.=FALSE)
}
if (is.null(opt$sample_ids)){
    print_help(opt_parser)
    stop("Please provide sample ids associated with peak files.", call.=FALSE)
}

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}

PeakFiles <- unlist(strsplit(opt$peak_files,","))
SampleIDs <- unlist(strsplit(opt$sample_ids,","))
if (length(PeakFiles) != length(SampleIDs)) {
    print_help(opt_parser)
    stop("Number of sample ids must equal number of peak files.", call.=FALSE)
}

## SCORE COLUMNS: data.frame with column, name and scale
columns <- data.frame(column=integer(), name=character(), scale=character())
for (spec in unlist(strsplit(opt$columns, ","))) {
    fields <- unlist(strsplit(spec, ":", fixed=TRUE))
    if (length(fields) < 2 || is.na(suppressWarnings(as.integer(fields[1])))) {
        stop(paste("Invalid --columns entry:", spec), call.=FALSE)
    }
    scale <- if (length(fields) > 2) fields[3] else ''
    if (!scale %in% c('', 'log2', 'log10', 'neglog10')) {
        stop(paste("Unknown scale in --columns entry:", spec), call.=FALSE)
    }
    columns <- rbind(columns, data.frame(column=as.integer(fields[1]), name=fields[2], scale=scale))
}

################################################
################################################
## READ IN DATA                               ##
################################################
################################################

## Long format: one row per peak and measure
plot.dat <- data.frame()
summary.dat <- data.frame()
for (idx in 1:length(PeakFiles)) {

    sampleid <- SampleIDs[idx]
    peaks <- read.table(PeakFiles[idx], sep="\t", header=FALSE, quote="", comment.char="")

    measures <- list()
    for (i in seq_len(nrow(columns))) {
        values <- as.numeric(peaks[,columns$column[i]])
        if (columns$scale[i] == 'neglog10') {
            values <- pmin(-log10(values), 1000)
        }
        measures[[columns$name[i]]] <- values
    }
    measures[['length']] <- peaks[,3] - peaks[,2]

    ## GET SUMMARY STATISTICS
    for (cname in names(measures)) {
        sdat <- summary(measures[[cname]])
        sdat["num_peaks"] <- nrow(peaks)
        sdat["measure"] <- cname
        sdat["sample"] <- sampleid
        sdat <- t(data.frame(x=matrix(sdat),row.names=names(sdat)))
        summary.dat <- rbind(summary.dat,sdat)
        plot.dat <- rbind(plot.dat, data.frame(name=sampleid, measure=cname, value=measures[[cname]]))
    }
}
plot.dat$name <- factor(plot.dat$name, levels=sort(unique(as.character(plot.dat$name))))

SummaryFile <- file.path(opt$outdir,paste(opt$outprefix,".summary.txt",sep=""))
write.table(summary.dat,file=SummaryFile,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

################################################
################################################
## PLOTS                                      ##
################################################
################################################

## Sample labels are file names without spaces, so wrap them after the dots instead
wrap_axis_labels <- function(labels, width = 25) {
    sapply(labels, function(label) {
        parts <- regmatches(label, gregexpr("[^.]+\\.?", label))[[1]]
        lines <- c()
        line <- ""
        for (part in parts) {
            if (nchar(line) > 0 && nchar(line) + nchar(part) > width) {
                lines <- c(lines, line)
                line <- part
            } else {
                line <- paste0(line, part)
            }
        }
        paste(c(lines, line), collapse="\n")
    }, USE.NAMES=FALSE)
}

plot.theme <- theme(legend.position="none",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.text.y = element_text(colour="black"),
                    axis.text.x= element_text(colour="black",face="bold"),
                    axis.line.x = element_line(linewidth = 1, colour = "black", linetype = "solid"),
                    axis.line.y = element_line(linewidth = 1, colour = "black", linetype = "solid"))

## RETURNS VIOLIN PLOT OBJECT
violin.plot <- function(plot.dat,ylab,title,log) {

    plot  <- ggplot(plot.dat, aes(x=name, y=value)) +
                geom_violin(aes(colour=name,fill=name), alpha = 0.3) +
                geom_boxplot(width=0.1) +
                scale_x_discrete(labels = wrap_axis_labels) +
                xlab("") +
                ylab(ylab) +
                ggtitle(title) +
                plot.theme
    if (log == 'log10') {
        plot <- plot + scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
    }
    if (log == 'log2') {
        plot <- plot + scale_y_continuous(trans='log2',breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))
    }
    return(plot)
}

############################

PlotFile <- file.path(opt$outdir,paste(opt$outprefix,".plots.pdf",sep=""))
pdf(PlotFile,height=6,width=max(6,3*length(unique(plot.dat$name))))

## PEAK COUNT PLOT
peak.count.dat <- data.frame(name=factor(SampleIDs, levels=levels(plot.dat$name)), count=as.vector(table(plot.dat$name[plot.dat$measure == 'length'])[SampleIDs]))
plot  <- ggplot(peak.count.dat, aes(x=name, y=count)) +
            geom_bar(stat="identity",aes(colour=name,fill=name), position = "dodge", width = 0.8, alpha = 0.3) +
            scale_x_discrete(labels = wrap_axis_labels) +
            xlab("") +
            ylab("Number of peaks") +
            ggtitle("Peak count") +
            plot.theme +
            geom_text(aes(label = count, x = name, y = count), position = position_dodge(width = 0.8), vjust = -0.6)
print(plot)

## VIOLIN PLOTS
print(violin.plot(plot.dat=plot.dat[plot.dat$measure == 'length',],ylab="Peak length (bp)",title="Peak length distribution",log='log10'))
for (i in seq_len(nrow(columns))) {
    log <- if (columns$scale[i] %in% c('log2', 'log10')) columns$scale[i] else ''
    print(violin.plot(plot.dat=plot.dat[plot.dat$measure == columns$name[i],],ylab=columns$name[i],title=paste(columns$name[i],"distribution"),log=log))
}
dev.off()

################################################
################################################
################################################
################################################
