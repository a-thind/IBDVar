# S03_ibd_filter.R - filter by ibd regions
# load library
library(vcfR)
library(dplyr)

cat("\nScript: S03_annotate_IBD_sv.R\n")
cat(date(), "\n")

# files
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript s03_ibd_filter.R in_vcf.vcf out_dir annot_vars ibd_regions")
}

if (!file.exists(args[1])) {
  stop("input VCF file path does not exist.")
} else if (substr(args[1], nchar(args[1]) - 3, nchar(args[1]))!=".vcf") {
  stop("input file is not a VCF file")
} else {
  in_vcf <- args[1]
}

if (!dir.exists(args[2])) {
  stop("output folder path does not exist.")
} else {
  out_dir <- args[2]
}

if (!file.exists(args[3])) {
  stop("input SV overlaps with CCDS BED file (sv_ccds.bed) does not exist.")
} else if (substr(args[3], nchar(args[3]) - 3, nchar(args[3]))!=".bed") {
  stop("input file is not a BED file")
} else {
  annot_vars <- args[3]
}

if (!file.exists(args[4])) {
  stop("input SV overlaps with IBD regions BED file (ibd_overlaps.bed) does not exist.")
} else if (substr(args[4], nchar(args[4]) - 3, nchar(args[4]))!=".bed") {
  stop("input file is not a BED file")
} else {
  ibd_regions <- args[4]
}

cat("\nInput VCF: ", in_vcf)
cat("\nOutput folder: ", out_dir)
cat("\nSV-CCDS overlaps BED file: ", annot_vars)
cat("\nSV-IBD region overlaps BED file: ", ibd_regions)

# read vcf file
vcf <- read.vcfR(in_vcf)
tidy_vcf <- vcfR2tidy(vcf)
variants <- tidy_vcf$fix

# remove genotypes for now


colnames <- c("CHROM", "START", "ID", "REF", "ALT", "QUAL", "FILTER","INFO", "FORMAT",
              "CDS_CHROM", "CDS_START", "CDS_END", "NC_ACCESSION","GENE", 
              "GENE_ID", 'CDS_ID', 'CCDS_STATUS', 'STRAND','CDS_LOCATIONS', 
              'MATCH_TYPE', "OVERLAP")

ccds_sv_df <- read.delim(annot_vars, header=F)
ibd_sv_df <- read.delim(ibd_regions, header=F)

ibd_ccds <- ccds_sv_df[,3] %in% ibd_sv_df[,3]
ibd_filtered_vars <- ccds_sv_df[ibd_ccds,]

# create mask for genotypes
not_genotypes <- sapply(ibd_filtered_vars, function(x){
  if (!any(grepl("/", x))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

# filter out genotypes
ibd_filtered_vars <- ibd_filtered_vars[,not_genotypes]
# add column names
colnames(ibd_filtered_vars) <- colnames

ibd_filtered_vars <- as_tibble(ibd_filtered_vars)

# split the info column by SVTYPE, 
filtered_vars <- ibd_filtered_vars %>% 
  mutate(SV_TYPE=gsub("(.*SVTYPE=)|(;.*)", "", INFO)) %>%
  mutate(SV_LENGTH=gsub("(.*SVLEN=)|(;.*)", "", INFO)) %>%
  mutate(SV_LENGTH=gsub(".*=.*", "", SV_LENGTH)) %>%
  mutate(END=gsub("(.*)(^END=)|(;.*)", "",INFO)) %>%
  mutate(END=gsub(".*=.*", "",END)) %>%
  mutate(CIGAR=gsub("(.*CIGAR=)|(;.*)", "", INFO)) %>%
  mutate(CIGAR=gsub(".*=.*", "", CIGAR)) %>%
  mutate(CI_POS=gsub("(.*CIPOS=)|(;.*)", "", INFO)) %>%
  mutate(CI_POS=gsub(".*=.*", "", CI_POS)) %>%
  mutate(CI_END=gsub("(.*CIEND=)|(;.*)", "", INFO)) %>%
  mutate(CI_END=gsub(".*=.*", "", CI_END))

# write filtered variants to a tab-delimited text file
write.table(filtered_vars, file=file.path(out_dir,"ibd_annotated_sv.tsv"), 
            sep="\t", row.names = F, quote=F)

cat("\nDone.\n")
