#!/bin/bash

# Anisha Thind, 29May2022

# Intended use:
# ./s01_add_variant_ids.sh $VCF &> s01_add_variant_ids.log

# Adds ID for each variant

# stop at runtime errors
set -e

# Starting message
echo $0
date
echo ""

# Progress report
bcftools --version
date
echo ""

# set files and folders
VCF=$1
outdir="/home/share/data/s04_annotate"
mkdir -p "${outdir}"
filepath=` basename "${VCF}" .vcf.gz `
basename="${filepath%%.*}"
out_vcf="${outdir}/${basename}".ID.vcf.gz

# input VCF
echo ""
echo "Input VCF: ${1}"
echo "Output VCF: ${out_vcf}"
echo ""

# Add variant IDs
echo "Adding variant IDs..."
bcftools annotate "${VCF}" \
  -I '%CHROM\_%POS\_%REF\_%ALT' \
  --threads 4 \
  -Oz -o "${out_vcf}"

# index VCF
echo "Indexing VCF file..."
bcftools index "${out_vcf}"
echo ""

# Completion message
echo "Done."
date
echo ""
