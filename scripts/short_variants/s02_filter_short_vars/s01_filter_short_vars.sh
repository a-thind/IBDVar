#!/bin/bash
# s01_filter_short_vars.sh - filter variants so only those that pass all filters and
#   meet DP and GQ thresholds are retained.
# Anisha Thind, 16May2022

# Example:
# ./s01_retain_pass_filter_vars.sh in_vcf out_dir &> s01_retain_pass_filter_vars.sh
# Parameters:
#   in_vcf: input VCF file
#   out_dir: output directory
#   GQ: geotype quality
#   DP: depth
#   threads: number of threads

# Stop at runtime errors, if any variable value is unset or a non-zero pipe status
set -euo pipefail

printf "Script: s01_filter_short_vars.sh\n"
date
echo ""

# VCF file input
in_vcf="${1}"
out_dir="${2}"
GQ="${3}"
DP="${4}"
threads="${5}"

# create directory for output
basename=$( basename "${in_vcf}" .vcf.gz ) 
out_dir="${out_dir}/s02_filter_short_vars"
stats_dir="${out_dir}/bcfstats"
mkdir -p "${out_dir}"
mkdir -p "${stats_dir}"

# name for filtered VCF
pass_vcf="${out_dir}/${basename}.pass_filtered.vcf.gz"
gq_dp_vcf="${pass_vcf%.pass_filtered*}.gq_dp.vcf.gz"
filtered_vcf="${out_dir}/${basename}.filtered.vcf.gz"

# check VCF file exists
if [ -z "${in_vcf}" ]; then
  echo "Missing argument: VCF file"
  exit 1
elif [ ! -e "${in_vcf}" ]; then
  echo "Error: Input VCF file ${in_vcf} not found."
  exit 1
fi

if [ -z "${out_dir}" ]; then
  echo "Error: Missing argument for output directory."
  exit 1
elif [ ! -e "${out_dir}" ]; then
  echo "Error: ${out_dir} is not a directory."
  exit 1
fi

# Retain variants that pass all filters
echo "Filtering variants, retaining only those that pass all filters..."
bcftools view -Oz -f "PASS" "${in_vcf}" > "${pass_vcf}"
echo ""

bcftools query -f "%FILTER\n" "${in_vcf}" | sort | uniq -c
echo ""

bcftools view -H "${in_vcf}" \
  | awk 'END{printf("Number of variants before filtering: %s\n", NR)}' 
bcftools view -H "${pass_vcf}" \
  | awk 'END{printf("Number of variants after filtering (PASS): %s\n\n", NR)}'

printf "Filtering variants with GQ>= %s and DP >= %s in all samples...\n" \
  "${GQ}" "${DP}"

bcftools view -i "GQ>=${GQ} && FORMAT/DP>=${DP}" "${pass_vcf}" \
  --threads "${threads}" \
  -Oz \
  -o "${filtered_vcf}" 

zgrep -v "^#" "${filtered_vcf}" \
  | awk 'END{printf("Number of variants after filtering: %s\n\n", NR)}' 

# Index VCF file
echo "Indexing Filtered VCF file..."
bcftools index "${filtered_vcf}"
echo ""

# completion message
echo "Done."
date
echo ""
