#!/bin/env Rscript
# s02_filter_ccds.R - filters CCDS data so transcripts with longest sequence for a
#   gene are retained only.
# Anisha Thind, 3Jul2022

# install packages if required
if (!require(BiocManager)) {
    install.packages("BiocManager")
}
if (!require(GenomicRanges)) {
    BiocManager::install("GenomicRanges")
}

# load libraries
library(GenomicRanges)

# clear workspace
rm(list=ls())
graphics.off()

file="ccds/CCDS.current.bed"
out_file="ccds/CCDS.current_filtered.bed"

# read ccds data
cds <- read.table(file, sep='\t', header=F)


cat(sprintf("Total number of coding sequences in CCDS file: %s\n", 
    nrow(cds)))

# create GRanges object
cdsRanges <- GRanges(seqnames=cds[,1], ranges=IRanges(start=cds[,2], 
    end=cds[,3]), strand=cds[,9], nc_accession=cds[,4], 
    gene=cds[,5], gene_id=cds[,6], cds_id=cds[,7], 
    ccds_status=cds[,8], cds_locations=cds[,10], 
    match_type=cds[,11])

# find longest transcripts for each gene
cdsMax <- aggregate(width(cdsRanges) ~ gene, data=cdsRanges, max)
# subset these maximum length rows
cds$width <- width(cdsRanges)
colnames(cdsMax)[2] <- "width"
filtered_cds <- merge(cdsRanges, cdsMax, by=c("gene", "width"))

# rearrange columns
filtered_cds <- filtered_cds[, c('seqnames', 'start', 'end', 'nc_accession', 
    'gene', 'gene_id', 'cds_id', 'ccds_status', 'strand','cds_locations', 
    'match_type')]

cat(sprintf("Total number of coding sequences in filtered CCDS file: %s\n", 
    nrow(filtered_cds)))
# write out filtered bed file
write.table(filtered_cds, file=out_file, quote=F, sep='\t', col.names=F )