#!/bin/bash

# Anisha Thind, 16May2022

# Intended use:
# ./s01_retain_pass_filter_vars.sh $VCF &> s01_retain_pass_filter_vars.sh
# Filters VCF retaining variants that pass all filters

# Stop at runtime errors
set -e

echo $0
date
echo ""

# VCF file input
VCF=$1

# create directory for output
basename=`basename "${VCF}" .vcf.gz`
data_dir="${VCF%%/`basename "${VCF}"`}"
out_dir="/home/share/data/s02_retain_pass_filter_vars"
mkdir -p "${out_dir}"
# name for filtered VCF
PASS_VCF="${out_dir}/${basename}.pass_filtered.vcf.gz"

# check VCF file exists
if [ -z "${VCF}" ]
then
  echo "Missing argument: VCF file"
  exit 1
fi
if [ ! -e "${VCF}" ]
then
  echo "Error: File ${VCF} cannot be found."
  echo "Terminating script..."
  exit 1
fi

# Retain variants that pass all filters
echo "Filtering variants, retaining only those that pass all filters..."
bcftools view -Oz -f "PASS" "${VCF}" > "${PASS_VCF}"
echo ""

# Index VCF file
echo "Indexing Filtered VCF file..."
bcftools index "${PASS_VCF}"
echo ""

# completion message
echo "Done."
date
echo ""
