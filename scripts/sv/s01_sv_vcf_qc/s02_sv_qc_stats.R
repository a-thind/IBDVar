#!/usr/bin/env Rscript
# s02_sv_qc_stats.R - computes summary stats for SV VCF

# load libraries
library(vcfR)
library(RColorBrewer)

source("utils/stats.R")

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

cat("\nInput VCF file:", in_vcf)
cat("\nOutput folder:", out_dir)

# read vcf file
vcf <- read.vcfR(in_vcf)
tidy_vcf <- vcfR2tidy(vcf)
meta <- tidy_vcf$meta
variants <- tidy_vcf$fix

# chromosome summary
chrom_plot <- variants %>% counts(CHROM) %>% 
  ggplot(aes(x=gsub("chr", "", CHROM), y=count, fill="")) + 
  geom_bar(stat="identity") +
  xlab("Chromosome") +
  ylab("Count") +
  theme(legend.position="none")
ggsave(file.path(out_dir, "chromosome_variants.png"), chrom_plot)
# sv types
variants %>% counts(SVTYPE)
type_plot<- variants %>% counts(SVTYPE) %>% plot_counts("Type")
ggsave(file.path(out_dir, "sv_type.png"), type_plot)

# Quality
qual_plot <- variants %>% filter(QUAL < 100) %>% ggplot(aes(x=QUAL)) +
  geom_histogram(bins=40, fill=rgb(0,0,1,0.5)) +
  scale_x_continuous(limits=c(0, 100)) +
  ggtitle("QUAL Score distribution < 100")
ggsave(file.path(out_dir, "sv_qual.png"), qual_plot)

cat("\nVariant Summary:\n")
summary(variants$QUAL)

# SV length
variants <- variants %>% mutate(across(SVLEN, as.numeric))
cat("SV length summary:\n")
summary(variants$SVLEN)
png(file.path(out_dir, "sv_len.png"))
hist(table(abs(variants$SVLEN)), xlim=c(0, 100), 
                    xlab="SV Length", col=rgb(0,0,1,0.5), 
                    main="SV Length Distribution < 100")
dev.off()
# imprecise
table(variants$IMPRECISE)

# filters
# all filters
filters_plot <- variants %>% counts(FILTER) %>% plot_counts(
  var_name = "Filter", flip=TRUE)
ggsave(file.path(out_dir,"filters.png"), filters_plot)
# pass filters by type
pass_plot <- variants %>% filter(FILTER=="PASS") %>% counts(SVTYPE) %>% 
  plot_counts(var_name = "Type", title="Variants passing all filters (PASS)")
ggsave(file.path(out_dir, "pass_filters.png"), pass_plot)
# pass out of total variants

cat("\nAnalysis Complete.\n")
