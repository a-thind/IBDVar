#!/bin/env Rscript
# s07_draw_ideogram.R - draws an ideogram from IBIS output segment file
# Anisha Thind, 30Jun2022 

# install packages

if (!require("optparse")) {
    install.packages("optparse")
}
if (!require("DescTools")) {
    install.packages("DescTools")
}


# load libraries
library(optparse)
library(scales)
library(ideogram)
library(DescTools)
library(RColorBrewer)

# clear workspace
rm(list=ls())
graphics.off()

# declare options for script
option.list <- list(
    make_option(
        c('-f', '--file', type='character', metavar='character', 
        default=NULL, help='Input IBIS IBD segments file')),
    make_option(
        c('-o', '--outdir', type='character', metavar='character', 
            default=NULL, help='Output directory'))
)

# read the command line args
opt.parser <- OptionParser(option_list=option.list)
# parse args
opt <- parse_args(opt.parser)

# check parsed options
if (is.null(opt$file)) {
    stop("Input file has not been provided.")
}

if (is.null(opt$outdir)) {
    stop("Output directory has not been provided.")
}

if (!file.exists(opt$file)) {
    stop(sprintf("File '%s' does not exist.", opt$file))
}

if (!file.exists(opt$outdir)) {
    stop(sprintf("Output folder '%s' does not exist.", opt$outdir))
}

ibd <- read.table(opt$file, header=TRUE)

# concatenate two sample names together
ibd$name <- paste0(ibd$sample1, '-', ibd$sample2)
# remove first two sample columns
ibd <- ibd[, -c(1,2)]

ibd$name <- as.factor(ibd$name)
# create colours
color <- hue_pal()(length(levels(ibd$name)))
# add alpha
color <- SetAlpha(color, alpha=0.65)
# map the colours to the samples
ibd$color <- factor(ibd$name, labels=color)

# create GRanges object


# subset required annotation columns
annots <- ibd[,c("name", "chrom", "phys_start_pos", "phys_end_pos", "color")]
colnames(annots)[2:4] <- c("chr","start", "stop")
annots$chr <- as.character(annots$chr)

# draw ideogram
p <- ideogram(annots)

# save plot
saveWidget(p, sprintf("%s/ideogram.html", opt$outdir), selfcontained=FALSE)
cat(sprintf("\nIdeogram HTML page (ideogram.html) created in '%s'.\n", opt$outdir))