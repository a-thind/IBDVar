# s07_draw_ideograms.R - draws ideograms from IBIS output segment file
# Anisha Thind, 30Jun2022 

# install packages
#devtools::install_github("freestatman/ideogRam")
install.packages("htmltools")

# load libraries
library(ideogRam)
library(GenomicRanges)

# clear workspace
rm(list=ls())
graphics.off()

file="/home/share/data/output/IHCAPX8/s07_select_haploblocks/ibis/ibis.seg"

ibd <- read.table(file, header=T)

# concatenate two sample names together
ibd$samples <- paste0(ibd$sample1, '-', ibd$sample2)
# remove first two sample columns
ibd <- ibd[, -c(1,2)]

ibd$samples <- as.factor(ibd$samples)
# create colours
colours <- hue_pal()(length(levels(ibd$samples)))
colours <- col2rgb(colours)
# map the colours to the samples
ibd$colours <- factor(ibd$samples, labels=colours)

# create GRanges object
ibd.seg <- GRanges(seqnames=ibd$chrom, 
    ranges=IRanges(start=ibd$phys_start_pos, end=ibd$phys_end_pos), 
    strand=rep('*', nrow(ibd)), 
    samples=ibd$samples,
    color=as.character(ibd$colours))

# draw ideogram
p <- ideogRam(organism="human") %>% 
    set_option(chrMargin=2, annotationsLayout = "overlay") %>% 
    set_option(showAnnotTooltip = TRUE) %>%
    add_track(ibd.seg)
# save plot
htmltools::save_html(p, "/home/share/data/output/ideogram.html")