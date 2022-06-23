# s08 explore IBD regions
# Anisha Thind, 22Jun2022

# Install packages
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.refGene")

# load libraries
library(GenomicRanges)

# clear workspace
rm(list=ls())
graphics.off()

# read IBD segment file from IBIS
ibis <- read.table("../../data/s07_select_haploblocks/ibis/IHCAPX8_ibis.seg", 
                   header=T)
# read Truffle segment file
truffle <- read.table("../../data/s07_select_haploblocks/truffle/IHCAPX8_vcf_truffle.segments", 
                      header=T)

ibisGr <- GRanges(seqnames=ibis$chrom, ranges=IRanges(start=ibis$phys_start_pos, 
                                                      end=ibis$phys_end_pos))
ibisGr
