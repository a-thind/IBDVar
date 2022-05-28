#!/bin/bash

# s02_check_vcf_stats.sh
# Anisha Thind, 17May2022

# Intended use:
# ./s02_check_vcf_stats.sh VCF_file &> s02_check_vcf_stats.log

# Checks the VCF stats of filtered VCF file (all PASS)

# stop at runtime errors
set -e

# start message
echo $0
date
echo ""

# set files and folder variables
VCF=$1
basename=` basename "${VCF}" .vcf.gz `
data_dir="/home/share/data"
stats_dir="${data_dir}/s02_retain_pass_filter_vars/bcfstats"
stats_file="${stats_dir}/${basename}.vchk"
mkdir -p "${stats_dir}"

echo "Check VCF stats for filtered VCF (all PASS filters)"
echo ""

# Generate progress report
bcftools --version
date
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
