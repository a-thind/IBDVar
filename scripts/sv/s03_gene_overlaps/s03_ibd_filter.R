# S03_ibd_filter.R - filter by ibd regions
# load library
library(vcfR)
library(tidyr)

rm(list=ls())

# files
base_dir="/run/user/1000/gvfs/sftp:host=138.250.31.2,user=anisha/"
in_vcf=file.path(base_dir, '/home/share/data/output/IHCAPX8/sv/s02_filter_sv/IHCAPX8_SV_dragen_joint.sv.pass.read_support.vcf')
annot_vars=file.path(base_dir, '/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/sv_ccds.bed')
ibd_regions=file.path(base_dir, '/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/ibd_overlaps.bed')
out_dir=file.path(base_dir, '/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/')
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
write.table(filtered_vars, file=file.path(out_dir,"ibd_annotated_sv.txt"), 
            sep="\t", row.names = F, quote=F)

