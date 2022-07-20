#!/bin/env Rscript
# Anisha Thind, 9Jul2022

# install packages
install.packages("rlang")
install.packages("vcfR")

# load libraries
library(vcfR)
library(dplyr)
library(ggplot2)



# clear workspace
rm(list=ls())
graphics.off()

getwd()

in_vcf = "/home/share/data/output/IHCAPX8/s05_select_variants/IHCAPX8_dragen_joint.clinvar.reheaded.VEP.AF.vcf.gz"

# read data into R
cat("Reading VCF into R...")
in_vcf <- file.path(in_vcf)
vcf <- read.vcfR(in_vcf, verbose=F)
vcf

cat("Converting to tidy format...")
vcf_tidy <- vcfR2tidy(vcf, info_only=TRUE, dot_is_NA=TRUE)
variants <- vcf_tidy$fix
variants
# extract metadata
meta <- vcf_tidy$meta
filtered_vars <- variants %>% filter(QUAL > 30)
filtered_vars
# plot histogram of QUAL
filtered_vars %>% filter(QUAL < 300) %>% ggplot(aes(x=QUAL)) + 
  geom_histogram(bins=40) + 
  scale_x_continuous(limits=c(0, 300), expand=c(0, 5)) + 
  ggtitle("QUAL score distribution < 300")

# extract benign assertions from clinvar
benign <- filtered_vars %>% select(CLNSIG, matches("benign|Benign"))
filtered_vars %>% 
  group_by(CLNREVSTAT) %>% 
  summarise(count=n()) %>% arrange(desc(count)) 
clinvar_benign <- filtered_vars %>% filter(grepl("benign|Benign", CLNSIG))
clinvar_benign_conf <- clinvar_benign %>% 
  group_by(CLNREVSTAT) %>% summarise(count=n())

# extract high impact assertions from VEP
vep_impact <- filtered_vars %>% group_by(vep_IMPACT) %>% 
  summarise(count=n()) %>% arrange(desc(count))

vep_impact_benign <- filtered_vars %>% filter(grepl("LOW", vep_IMPACT))
vep_impact

# cadd below 20
cadd <- filtered_vars %>% filter(vep_CADD_PHRED >= 20)
cadd_benign <- filtered_vars %>% filter(vep_CADD_PHRED < 20)
cadd_benign
# filter sift
polyphen_benign <- filtered_vars %>% filter(grepl("benign", vep_PolyPhen))


sift <- filtered_vars %>% select(vep_SIFT)%>% distinct()
sift_benign <- filtered_vars %>% filter(grepl("tolerated", vep_SIFT))
sift_benign

filtered_vars %>% group_by(CLNREVSTAT) %>% summarise()

vep_consequnce <- filtered_vars %>% distinct(vep_Consequence)


test_vars <- filtered_vars %>% 
  filter(!grepl("benign|Benign", CLNSIG)) %>%
  filter(vep_CADD_PHRED >= 20) %>%
  filter(!grepl("tolerated", vep_SIFT)) %>%
  filter(!grepl("benign", vep_PolyPhen)) %>%
  filter(!grepl("LOW", vep_IMPACT)) %>%
  filter(!grepl("synonymous", vep_Consequence))

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

filtered_vars %>% counts(CLNSIG) 
filtered_vars %>% counts(CLNSIG) %>% plot_counts()
filtered_vars %>% counts(vep_SIFT)
filtered_vars %>% filter(grepl("[a-z]+", vep_SIFT))
