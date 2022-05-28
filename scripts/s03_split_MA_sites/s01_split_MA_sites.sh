#!/bin/bash

# s01_split_MA_sites.sh
# Anisha Thind, 18May2022

# Intended use:
# ./s01_split_MA_sites.sh filtered_vcf &> s01_split_MA_sites.log
# Split alternate multiallele sites

# stop at runtime errors
set -e

# start message
echo $0
date
echo ""

# set file and folder variables
VCF=$1
basename=` basename ${VCF} .vcf.gz `
data_dir="/home/share/data"
out_dir="${data_dir}/s03_split_MA_sites"
mkdir -p "${out_dir}"
MA_VCF="${out_dir}/${basename}.split_MA.vcf.gz"

# Check VCF argument 
if [ ! -z "${VCF}" ]
then
    echo "Error: missing VCF file argument."
    echo "Aborting script..."
    exit 1
fi

# Make progress report
bcftools --version
date
echo ""

# Split MA sites
echo "Split multiallelic sites..."
bcftools norm -m-any "${VCF}" -Oz -o "${MA_VCF}"

# Indexing
echo "Indexing VCF..."
bcftools index "${MA_VCF}"
echo ""

# Adding BCF stats counts

# Completion message
echo "Done."
date
echo ""
