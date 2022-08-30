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

colnames <- c("CHROM", "START", "ID", "REF", "ALT", "QUAL", "FILTER","INFO", "FORMAT",
              "CDS_CHROM", "CDS_START", "CDS_END", "NC_ACCESSION","GENE", 
              "GENE_ID", 'CDS_ID', 'CCDS_STATUS', 'STRAND','CDS_LOCATIONS', 
              'MATCH_TYPE', "OVERLAP")

ccds_sv_df <- read.delim(annot_vars, header=F)
ibd_sv_df <- read.delim(ibd_regions, header=F)

ibd_ccds <- ccds_sv_df[,3] %in% ibd_sv_df[,3]
ibd_filtered_vars <- ccds_sv_df[ibd_ccds,]

head(ibd_filtered_vars)

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

# aggregate variants by ID and get a list of overlapping genes per variant
grouped <- ibd_filtered_vars %>% group_by(ID) %>% summarise(GENES=toString(GENE))

ibd_filtered_vars <- variants %>% filter(ID %in% grouped$ID)

ibd_annot_vars <- left_join(ibd_filtered_vars, grouped, by=c("ID"="ID"))

# rename position column as start
ibd_annot_vars <- ibd_annot_vars %>% rename(START=POS) %>% 
                                      mutate(SVLEN=abs(SVLEN))



# write filtered variants to a tab-delimited text file
write.table(ibd_annot_vars, file=file.path(out_dir,"ibd_annotated_sv.tsv"), 
            sep="\t", row.names = F, quote=F)

cat("\nDone.\n")
