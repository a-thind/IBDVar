#!/bin/bash
# s02_check_vcf_stats.sh - Checks the VCF stats of filtered VCF file (all PASS)
# Anisha Thind, 17May2022

# Intended use:
# ./s02_check_vcf_stats.sh out_dir &> s02_check_vcf_stats.log
# Parameters:
#   out_dir: output directory

# stop at runtime errors
set -e
# stop if any variable value is not set
set -u
# stop pipeline if non-zero status
set -o pipefail

# start message
echo -e "Script: s02_check_vcf_stats.sh"
date
echo ""

# set files and folder variables
out_dir="${1}/s02_retain_pass_filter_vars"
filtered_vcf=$( find "${out_dir}" -name *.filtered.vcf.gz ) 

if [ ! -e "${filtered_vcf}" ]; then
  echo "Error: filtered VCF file not found."
  exit 1
fi

basename=$( basename "${filtered_vcf}" .pass_filtered.vcf.gz ) 
stats_dir="${out_dir}/bcfstats"
stats_file="${stats_dir}/${basename}.vchk"
mkdir -p "${stats_dir}"

echo -e "Checking VCF stats for filtered VCF...\n"

# Generate progress report
bcftools --version
date
echo ""

echo "Names of samples:"
bcftools query -l "${filtered_vcf}"
echo ""

# count variants
echo "Variant Counts"
echo "-----------------"
bcftools +counts "${filtered_vcf}"
echo ""

bcftools view -H "${filtered_vcf}" \
  | awk 'END{ printf("Number of variants retained (PASS): %d\n\n", NR) }'

# generate BCFstats
echo "Calculating BCF stats..."
bcftools stats -s - "${filtered_vcf}" > "${stats_file}"
echo ""

# generate BCF stats plot
echo "Generating BCF stats plots..."
plot-vcfstats -p "${stats_dir}" "${stats_file}"
echo ""

# Completion message
echo "Done."
date
echo ""
