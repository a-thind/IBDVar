# S03_ibd_filter.R - filter by ibd regions
# load library
library(vcfR)
library(GenomicRanges)
library(rtracklayer)

rm(list=ls())

in_vcf='/run/user/1000/gvfs/sftp:host=138.250.31.2,user=anisha/home/share/data/output/IHCAPX8/sv/s02_filter_sv/IHCAPX8_SV_dragen_joint.sv.pass.read_support.vcf.gz'
ibd_bed='/run/user/1000/gvfs/sftp:host=138.250.31.2,user=anisha/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/ibd.bed' 
out_vcf='/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/ibd_filtered.vcf'

cat("\nReading in VCF...")
vcf <- read.vcfR(in_vcf)
tidy_vcf <- vcfR2tidy(vcf)
tidy_vcf

tidy_vcf$meta %>% print(n=30)


variants <- tidy_vcf$fix
ibd_df <- read.delim(ibd_bed, sep="\t", col.names=c("chrom", "start", "end"))

# if variant start position and end position in ibd regions
ibd_ranges <- GRanges(seqname=ibd_df$chrom, 
                      ranges=IRanges(start=ibd_df$start,
                                     end = ibd_df$end), 
                      strand=rep("*", nrow(ibd_df)))
ibd_ranges

# extract break ends and non break ends separately 
# break ends (do not have an end position)
bnd <- variants %>% filter(SVTYPE=="BND")
not_bnd <- variants %>% filter(SVTYPE!="BND")
# not bnd ranges
not_bnd_ranges <- GRanges(seqnames = not_bnd$CHROM, 
                      ranges=IRanges(start=not_bnd$POS, end=not_bnd$END,
                                     names=not_bnd$ID),
                      strand=rep("*", nrow(not_bnd)))
# find overlaps 
not_bnd_overlaps <- findOverlaps(not_bnd_ranges, ibd_ranges)
not_bnd_overlaps
# extract variant IDs of overlaps
ibd_not_bnd <- names(not_bnd_ranges)[queryHits(not_bnd_overlaps)]
# drop all non-bnd variants that are no overlaps with IBD regions
not_bnd_ranges <- not_bnd_ranges[names(not_bnd_ranges) %in% ibd_not_bnd]
cat("Number of non break-ends in IBD regions:", length(ibd_not_bnd))

# for break-ends find overlaps
bnd_ranges <- GRanges(seqname=bnd$CHROM, 
                      ranges=IRanges(start=bnd$POS, width=1, names=bnd$ID), 
                      strand=rep("*", nrow(bnd)))
bnd_overlaps <- findOverlaps(bnd_ranges, ibd_ranges)
# extract variant IDS of break-ends that overlap with IBD variants
ibd_bnd <- names(bnd_ranges)[queryHits(bnd_overlaps)]
# drop all bnds that are not present in IBD regions
bnd_ranges <- bnd_ranges[names(bnd_ranges) %in% ibd_bnd,]

cat("Number of break-ends in IBD regions:", length(ibd_bnd))
# all variants in IBD ranges
filtered_vars <- variants %>% filter(ID %in% c(ibd_not_bnd, ibd_bnd))
# read in CCDS bed
ccds <- vroom::vroom('/home/share/resources/ccds/ccds_filtered.bed', 
                     delim="\t", col_names = c('seqnames', 'start', 'end', 'nc_accession', 
                    'gene', 'gene_id', 'cds_id', 'ccds_status', 'strand',
                    'cds_locations','match_type'))
# create ranges object for CCDS
ccds_ranges <- GRanges(seqnames=ccds$seqnames, 
                       ranges = IRanges(start=ccds$start, 
                                        end=ccds$end,
                                        names=ccds$cds_id),
                       strand=ccds$strand)

bnd_ccds_overlaps <- findOverlaps(bnd_ranges, ccds_ranges)
# break-ends IDs with ccds overlaps
bnd_ccds <- names(bnd_ranges)[queryHits(bnd_ccds_overlaps)]
# non-bnd IDs with ccds overlaps
not_bnd_ccds_overlaps <- findOverlaps(not_bnd_ranges, ccds_ranges)
not_bnd_ccds <- names(not_bnd_ranges)[queryHits(not_bnd_ccds_overlaps)]
not_bnd_ccds_overlaps
bnd_ccds_overlaps

# get break end overlap cds ids
bnd_ccds_cds <- names(ccds_ranges)[subjectHits(bnd_ccds_overlaps)]
not_bnd_ccds_cds <- names(ccds_ranges)[subjectHits(not_bnd_ccds_overlaps)]
not_bnd_ccds
