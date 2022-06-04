#!/bin/bash

# s02_check_clinvar.vcf.sh

# Anisha Thind, 4Jun2022

# Intended use:
# s02_check_clinvar_vcf.sh $VCF &> s02_check_clinvar_vcf.log

# Stop at runtime errors
set -e

# Starting message
echo $0
date
echo ""

# files and folders
VCF=$1
clinvar="/home/share/resources/clinvar/clinvar_20220507.vcf.gz"

# contig names from clinvar vcf
bcftools --version
date
echo ""

# generate progress report
echo "Input VCF: ${1}"
echo "ClinVar VCF: ${clinvar}"
echo ""

# Index clinvar vcf
echo "Indexing ClinVar VCF..."
bcftools index -f "${clinvar}"

echo "Checking ClinVar VCF..."
# for each vcf get contig names
bcftools index "${clinvar}" -s | awk '{print $1}' | head -n 25 > clinvar_chr.txt
bcftools index "${VCF}" -s | awk '{print $1}' | head -n 25 > vcf_chr.txt

echo "Creating translation table for chromosomes..."
# combine chromosome names for translation table text file
paste -d"\t" clinvar_chr.txt vcf_chr.txt > s02_clinvar_chr_map.txt
echo "Chromosome translation table created."
# clean up contig names files
rm clinvar_chr.txt vcf_chr.txt
echo ""

echo "Reference in VCF data:"
echo ""
bcftools view -h ${VCF} | grep "^##reference"
echo ""

echo "Reference in ClinVar:"
echo ""
bcftools view -h "${clinvar}" | grep "^##reference"
echo ""

# Completion message
echo "Done."
date
echo ""

