#!/bin/bash

# s02_bcfstats.sh
# Anisha Thind, 10May2022

# Intended use:
# ./s02_bcfstats.sh vcf_file &> s02_bcfstats.log
# Performs QC for VCF files

# stop at runtime errors
set -e

# start message
echo $0
date
echo ""

# set variables for folders
VCF="$1"
data_dir="/home/share/data/s01_short_vcf_check/"
stats_dir="${data_dir}/bcfstats"

# check VCF file exists
if [ ! -e "${VCF}" ]
then
   echo "VCF file: ${VCF} not found."
   exit 1
fi

# create output files
base_name=` basename ${VCF} `
# make stats directory
mkdir -p "${stats_dir}"
stats_file="${stats_dir}/${base_name}.vchk"

# Generate progress report
bcftools --version
echo ""
echo "VCF file: ${VCF}"
echo "Data folder: ${data_dir}"
echo ""

# reindex vcf file
echo "Indexing VCF file..."
bcftools index -f "${VCF}"
echo ""

# count variants
echo "Variant Counts"
echo "-----------------"
bcftools +counts "${VCF}"
echo ""

# generate BCFstats
echo "Calculating BCF stats..."
bcftools stats -s - "${VCF}" > "${stats_file}"
echo ""

# generate BCF stats plot
echo "Generating BCF stats plots..."
plot-vcfstats -p "${stats_dir}" "${stats_file}"
echo ""

# Completion message
echo "Done."
date
echo ""
