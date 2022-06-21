#!/bin/bash
# s01_retain_pass_filter_vars.sh - uses "FILTER" field to retain only those variants
# which pass all filters (PASS)
# Anisha Thind, 16May2022

# Intended use:
# ./s01_retain_pass_filter_vars.sh in_vcf out_dir &> s01_retain_pass_filter_vars.sh
# in_vcf: input VCF file
# out_dir: output directory

# Stop at runtime errors
set -e
# stop if any variable value is unset
set -u
# stop pipeline if non-zero status
set -o pipefail

echo $0
date
echo ""

# VCF file input
in_vcf=$1
out_dir=$2

# create directory for output
basename=`basename "${in_vcf}" .vcf.gz`
out_dir="${out_dir}/s02_retain_pass_filter_vars"
mkdir -p "${out_dir}"
# name for filtered VCF
pass_vcf="${out_dir}/${basename}.pass_filtered.vcf.gz"

# check VCF file exists
if [ -z "${in_vcf}" ]
then
  echo "Missing argument: VCF file"
  exit 1
fi
if [ ! -e "${in_vcf}" ]
then
  echo "Error: File ${in_vcf} cannot be found."
  echo "Terminating script..."
  exit 1
fi

# Retain variants that pass all filters
echo "Filtering variants, retaining only those that pass all filters..."
bcftools view -Oz -f "PASS" "${in_vcf}" > "${pass_vcf}"
echo ""

# Index VCF file
echo "Indexing Filtered VCF file..."
bcftools index "${pass_vcf}"
echo ""

# completion message
echo "Done."
date
echo ""
