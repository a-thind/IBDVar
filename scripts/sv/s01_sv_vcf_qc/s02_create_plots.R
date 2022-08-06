#!/usr/bin/env Rscript

# load library
library(ggplot2)

setwd("/home/share/scripts/sv/s01_sv_vcf_qc")

file=

# read in data from CSV
df <- read.csv(file, header=F)

# remove surplus column
df <- df[-1]

# rename columns
colnames(df) <- c("Sample", "Type", "Number")

# remove extraneous details (e.g. "Number of" "(PASS)")
df[,2] <- sapply(df[,2], function(i) gsub("Number of | \\(PASS\\)", "", i))

# plot SV types by sample
ggplot(df, aes(fill=Type, y=Number, x=Sample)) + geom_bar(position="stack", 
                                                          stat="identity")

args = commandArgs(trailingOnly=TRUE)
# save plot
ggsave(filename="../../../data/sv/s01_sv_vcf_qc/sv_type_by_sample.png", 
       device="png")

# plot by SV type
ggplot(df, aes(fill=Sample, y=Number, x=Type)) + geom_bar(position="stack", 
                                                          stat="identity")

ggsave(filename="../../../data/sv/s01_sv_vcf_qc/sv_type_by_sample.png", 
       device="png")

