#!/bin/env Rscript
# s07_draw_ideogram.R - draws ideograms from IBIS output segment file
# Anisha Thind, 30Jun2022 

# install packages
#devtools::install_github("freestatman/ideogRam")
#install.packages("htmltools")
#install.packages("optparse")

# load libraries
library(ideogRam)
library(GenomicRanges)
library(optparse)
library(scales)

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

ibd <- read.table(opt$file, header=T)

# concatenate two sample names together
ibd$samples <- paste0(ibd$sample1, '-', ibd$sample2)
# remove first two sample columns
ibd <- ibd[, -c(1,2)]

ibd$samples <- as.factor(ibd$samples)
# create colours
colours <- hue_pal()(length(levels(ibd$samples)))
# map the colours to the samples
ibd$colours <- factor(ibd$samples, labels=colours)

# create GRanges object
ibd.seg <- GRanges(seqnames=ibd$chrom, 
    ranges=IRanges(start=ibd$phys_start_pos, end=ibd$phys_end_pos), 
    strand=rep('*', nrow(ibd)), 
    samples=ibd$samples,
    color=as.character(ibd$colours))

# draw ideogram
p <- ideogRam(organism="human", assembly="GRCh38") %>% 
    set_option(chrMargin=2, annotationsLayout = "overlay") %>% 
    set_option(showAnnotTooltip = TRUE) %>%
    add_track(ibd.seg)
# save plot
htmltools::save_html(p, sprintf("%s/ideogram.html", opt$outdir))
cat(sprintf("\nIdeogram HTML page (ideogram.html) created in '%s'.\n", opt$outdir))