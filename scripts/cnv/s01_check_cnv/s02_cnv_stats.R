#!/bin/env Rscript
# s01_cnv_stats.R

if (!require("vcfR") {
    install.packages("vcfR")
}

# load libraries
library(vcfR)

# clear environment
rm(list=ls())
graphics.off()

vcf_file <- "/home" 
vcf <- read.vcfR(vcf_file)