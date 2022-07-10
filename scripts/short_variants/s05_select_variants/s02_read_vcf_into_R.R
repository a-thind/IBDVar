#!/bin/env Rscript
# Anisha Thind, 9Jul2022

# load libraries
library(vcfR)

in_vcf = "../../../data/output/IHCAPX8/s05_filter_vep_vars/IHCAPX8_dragen_joint.clinvar.reheaded.VEP.AF.vcf.gz"

# read data into R
cat("Reading VCF into R...s")
vcf <- read.vcfR(in_vcf, verbose=F)
vcf

cat("Converting to tidy format...")
vcf_tidy <- vcfR2tidy(vcf, info_only=TRUE, dot_is_NA=TRUE)
vcf <- as.data.frame(meta$fix)