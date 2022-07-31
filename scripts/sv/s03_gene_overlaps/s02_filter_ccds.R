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
suppressMessages(suppressWarnings(library(GenomicRanges, quietly=TRUE)))
suppressMessages(suppressWarnings(library(optparse, quietly=TRUE)))

# clear workspace
rm(list=ls())
graphics.off()

opt_list <- list(
    make_option(c('-f', '--file', metavar='character', type='character',
        help='Input SV VCF file', default=NULL)),
    make_option(c('-o', '--outdir', metavar='character', type='character',
        help='Output folder for output overlaps file', default=NULL))
    )

# parse options
opt_parser <- OptionParser(option_list=opt_list)
opts <- parse_args(opt_parser)

# check parsed options
if (is.null(opts$file)) {
    stop("Input file has not been provided.")
}

if (is.null(opts$outdir)) {
    stop("Output directory has not been provided.")
}

if (!file.exists(opts$file)) {
    stop(sprintf("File '%s' does not exist.", opts$file))
}

if (!file.exists(opts$outdir)) {
    stop(sprintf("Output folder '%s' does not exist.", opts$outdir))
}

in_file <-  opts$file
out_file <- "ccds_filtered.bed"
out_file <- sprintf("%s/%s", opts$outdir, out_file)

cds <- read.table(in_file, sep='\t', header=F)


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
filtered_cds <- filtered_cds[!duplicated(filtered_cds$gene),]
cat(sprintf("Total number of coding sequences in filtered CCDS file: %s\n", 
    nrow(filtered_cds)))

# write out filtered bed file
write.table(filtered_cds, file=out_file, quote=F, sep='\t', col.names=F, row.names=F)
