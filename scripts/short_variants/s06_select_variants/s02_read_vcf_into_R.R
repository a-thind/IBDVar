#!/bin/env Rscript
# Anisha Thind, 9Jul2022

# install packages
if (!require("vcfR")) {
  install.packages("vcfR")
}
if (!require("ggvenn")) {
  install.packages("ggvenn")
}


# load libraries
library(vcfR)
library(dplyr)
library(ggplot2)
library(ggvenn)



# clear workspace
rm(list=ls())
graphics.off()

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error: both input VCF file and output folder arguments are required.")
}

if (!file.exists(args[1])) {
  stop("input VCF file path does not exist.")
} else if (substr(args[1], nchar(args[1]) - 6, nchar(args[1]))!=".vcf.gz") {
  stop("input file is not a VCF file")
} else {
  in_vcf <- args[1]
}

if (!dir.exists(args[2])) {
  stop("output folder path does not exist.")
} else {
  out_dir=args[2]
}

# output file name
out_file <- file.path(out_dir, "filtered_short_vars.tsv")

# starting message:
cat("\nSelecting protein-affecting variants\n")
cat(sprintf("Input VCF: %s\n", in_vcf))
cat(sprintf("Selected variants output file: %s\n\n", out_file))


#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
# Counts per feature function
# Example: variants_df %>% counts(feature)
# data: dataframe / tibble of variants
# feature: feature to aggregate by
counts <- function(data, feature) {
  data %>%
    group_by({{ feature }}) %>%
    summarise(count = n()) %>% arrange(desc(count))
}

# produce barplot showing counts for each type of group variable
plot_counts <- function(data) {
  ggplot(data, aes(x=pull(data, 1), y=count)) +
    geom_bar(stat="identity") +
    xlab("Type") +
    ylab("Count") +
    coord_flip()
}


# read data into R
cat("\nReading VCF into R...\n")
in_vcf <- file.path(in_vcf)
vcf <- read.vcfR(in_vcf, verbose=F)
vcf

cat("\nConverting to tidy format...\n")
vcf_tidy <- vcfR2tidy(vcf, info_only=TRUE, dot_is_NA=TRUE, format_types=TRUE)
variants <- vcf_tidy$fix
variants

# extract metadata
meta <- vcf_tidy$meta
meta %>% print(n=100)

# split polyphen and sift columns into score and call
variants <- variants %>%
  mutate(vep_SIFT=gsub(",\\.|\\.,", "", vep_SIFT)) %>%
  mutate(vep_SIFT_score=gsub("[a-z]+\\(|\\)","", vep_SIFT)) %>%
  mutate(vep_SIFT_call=gsub("[^a-z]+\\)", "", vep_SIFT)) %>%
  mutate(vep_SIFT_score=gsub("[a-z_]+", "", vep_SIFT_score)) %>%
  select(-vep_SIFT)

colnames(variants)

variants <- variants %>%
  mutate(vep_PolyPhen=gsub(",\\.|\\.,", "", vep_PolyPhen)) %>%
  mutate(vep_PolyPhen_score=gsub("[a-z]+\\(|\\)","", vep_PolyPhen)) %>%
  mutate(vep_PolyPhen_call=gsub("[^a-z]+\\)", "", vep_PolyPhen)) %>%
  mutate(vep_PolyPhen_score=gsub("[a-z_]+", "", vep_PolyPhen_score)) %>%
  select(-vep_PolyPhen)

# verify splitting plyphen and sift columns
variants %>% select(vep_SIFT_score, vep_SIFT_call, vep_PolyPhen_score,
                    vep_PolyPhen_call) %>%
  filter(vep_SIFT_score!=".") %>%
  head()

# summarise counts
variants %>% counts(vep_Consequence)
variants %>% counts(vep_SIFT_call)
variants %>% counts(vep_PolyPhen_call)
variants %>% counts(CLNSIG)
variants %>% filter(vep_CADD_PHRED >= 20) %>% 
  summarise(count=n())
variants %>% counts(vep_IMPACT)

# plot histogram of QUAL
variants %>% filter(QUAL < 100) %>% ggplot(aes(x=QUAL)) +
  geom_histogram(bins=40, fill=rgb(0,0,1,0.5)) +
  scale_x_continuous(limits=c(0, 100), expand=c(0, 5)) +
  ggtitle("QUAL Score distribution < 100")

# replace all dots with NAs
filtered_vars <- as.data.frame(variants)
filtered_vars[filtered_vars=="."] <- NA

# change cadd columns to numeric
filtered_vars$vep_CADD_PHRED <- as.numeric(filtered_vars$vep_CADD_PHRED)
filtered_vars$vep_CADD_RAW <- as.numeric(filtered_vars$vep_CADD_RAW)

# apply filters
# clinvar
clinvar <- filtered_vars$CLNSIG %in% c("Pathogenic/Likely_pathogenic")
# sift
sift <- filtered_vars$vep_SIFT_call %in% c("deleterious")
# PolyPhen
polyphen <- filtered_vars$vep_PolyPhen_call %in% c("probably_damaging")

# CADD > 20
cadd <- filtered_vars$vep_CADD_PHRED >= 20 &
  !is.na(filtered_vars$vep_CADD_PHRED)

summary(cadd)

impact <- filtered_vars$vep_IMPACT == "HIGH"

combined_filter <- clinvar |
  (sift & cadd | cadd & polyphen) | impact
cat(sprintf("Number of variants passing combined filter: %s\n\n", 
            sum(combined_filter)))

cat(sprintf("Number of LOF variants: %s\n", sum(impact)))

cat(
  sprintf("Number of variants identified as probably damaging (PolyPhen): %s\n\n", 
          sum(polyphen)))

cat(
  sprintf("Number of variants identified as deleterious (SIFT): %s\n\n", 
          sum(sift)))

cat(
  sprintf("Number of variants identified as pathogenic / likely pathogenic (ClinVar): %s\n\n", 
          sum(clinvar)))

cat(
  sprintf("Number of variants identified as protein-affecting (SIFT, PolyPhen, CADD): %s\n\n", 
          sum(sift & polyphen & cadd)))
cat(
  sprintf("Number of variants identified as deleterious (SIFT) & CADD Phred >= 20: %s\n\n", 
          sum(cadd & sift)))
cat(
  sprintf("Number of variants identified as probably damaging (PolyPhen) & CADD Phred >= 20: %s\n\n", 
          sum(cadd & polyphen)))

cat("\n")

hist(table(filtered_vars$vep_CADD_PHRED),
     main="CADD Phred Score Distribution",
     xlab = "CADD Phred Score", col=rgb(0,0,1,0.5))

# venn diagram between SIFT, PolyPhen and CADD
ggvenn(tibble("CADD"=cadd, "PolyPhen"=polyphen, "SIFT"=sift))
# venn diagram between SIFT, PolyPhen and CADD and VEP impact
ggvenn(tibble("CADD"=cadd, "PolyPhen"=polyphen, "SIFT"=sift, "LOF"=impact))
# venn diagram between SIFT, PolyPhen and CADD and VEP impact
ggvenn(tibble("CADD"=cadd, "PolyPhen"=polyphen, "SIFT"=sift, "clinvar"=clinvar))


all_filters_vars <- filtered_vars[combined_filter,]

# clean up consequence names
vars_table <- all_filters_vars %>%
  mutate(vep_Consequence=gsub("_variant", "", vep_Consequence)) %>%
  mutate(vep_Consequence=gsub("_", " ", vep_Consequence)) %>%
  mutate(vep_Consequence=gsub("&", ", ", vep_Consequence)) %>%
  mutate(vep_HGVSc=gsub(".*:", "", vep_HGVSc)) %>%
  mutate(vep_HGVSp=gsub(".*:", "", vep_HGVSp)) %>%
  mutate(vep_HGNC_ID=gsub(".*:", "", vep_HGNC_ID)) %>%
  mutate(vep_CLIN_SIG=gsub("&", "", vep_CLIN_SIG)) %>%
  mutate(CLNSIG=gsub("_", " ", CLNSIG)) %>%
  mutate(vep_SIFT_call=gsub("_", " ", vep_SIFT_call)) %>%
  mutate(vep_PolyPhen_call=gsub("_", " ", vep_PolyPhen_call)) %>%
  mutate(ID=sub("_", ":", ID)) %>%
  mutate(ID=gsub("_", " ", ID)) %>%
  # clean column names
  rename_with(~ toupper(gsub("vep_", "", .x, fixed = TRUE)))

write.table(vars_table, file=out_file, quote = F, row.names=F, sep="\t")
cat(sprintf("Number of protein-coding variants selected: %s\n\n", 
            nrow(vars_table)))

