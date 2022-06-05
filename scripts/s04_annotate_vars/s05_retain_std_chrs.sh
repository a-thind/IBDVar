#!/bin/bash

# Anisha Thind, 4Jun2022

# Intended use:
# ./s05_retain_std_chrs.sh $VCF &> s05_retain_std_chrs.sh

# stop at runtime errors
set -e

# start message
echo $0
date
echo ""

# files and folders
VCF=$1
std_chr="${VCF%%.vcf*}.std_chr.vcf.gz"
updated_vcf="${std_chr%%.std*}.reheaded.vcf.gz"
data_folder="/home/share/data/s04_annotate"
HEADER="${data_folder}/header.txt"

# progress report
bcftools --version
date
echo ""

echo "Input VCF file: ${VCF}"
echo "Output VCF file: ${updated_vcf}"
echo ""

echo "Selecting variants from standard chromosomes only..."
bcftools annotate "${VCF}" \
   -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
   --threads 4 \
   -Oz -o "${std_chr}"

# Index VCF
bcftools index "${std_chr}"
echo ""

# variant counts in original
echo "-- Variant counts in annotated VCF file --"
bcftools +counts "${VCF}"
echo ""
echo "-- Variant counts in standard chromosomes of annotated VCF --"
bcftools +counts "${std_chr}"
echo ""

# summary of variants per contig
echo "-- Variant counts in standard chromosomes --"
printf "%s\t%s\t\t%s\n" Contig Length Records
echo "------------------------------------"
bcftools index "${std_chr}" -s
echo ""

# update header of VCF
echo "Updating VCF header..."
bcftools view -h "${std_chr}"| zgrep -E '(^##contig=<ID=chr([0-9]{1,2}|X|Y),)|^##[^contig].*|^#CHROM.*' > "${HEADER}"
# reheader VCF 
bcftools reheader "${std_chr}" \
   -h "${HEADER}" \
   --threads 4 \
   -o "${updated_vcf}"
echo ""

# index VCF
echo "Indexing VCF file..."
bcftools index "${updated_vcf}"

# verify changes have been made
echo "Number of contig fields in annotated VCF file:"
bcftools view -h "${std_chr}" | grep "^##contig" | wc -l
echo ""
echo "Number of contig fields in standard chromosome VCF file:"
bcftools view -h "${updated_vcf}" | grep "^##contig" | wc -l
echo ""

# Completion message
echo "Done."
date
echo ""
