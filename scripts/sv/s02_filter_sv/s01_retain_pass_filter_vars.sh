#!/bin/bash
# s01_retain_pass_filter_var.sh - retain only variants that have passed all filters
# Anisha Thind, 25Jun2022
# Parameters:
#   $1: (in_vcf) input vcf
#   $2: (out_dir) output folder

# stop at runtime errors, any variable value is unset or non-zero status in pipe
set -euo pipefail

# files and folders
in_vcf="${1}"
out_dir="${2}/s02_filter_sv"
mkdir -p "${out_dir}"
basename=$( basename "${in_vcf}" .vcf.gz )
pass_vcf="${out_dir}/${basename}.pass.vcf.gz"

echo -e "Script: s01_retain_pass_filter_var.sh\n"
date
echo ""

# Progress report
bcftools --version
date
echo ""

echo "Input VCF: ${in_vcf}"
echo "PASS filtered VCF: ${pass_vcf}"
echo -e "Output folder: ${out_dir}\n"


bcftools view -H "${in_vcf}" | awk 'END{ printf("Total number of variants BEFORE filtering: %s\n\n", NR) }'  

echo -e "Filtering out variants that have not passed all filters..."
bcftools view -Oz -f "PASS" "${in_vcf}" > "${pass_vcf}"
echo ""

echo -e "Indexing VCF file...\n"
bcftools index "${pass_vcf}"

bcftools view -H "${pass_vcf}" | awk 'END{ printf("Total number of variants AFTER filtering: %s\n\n", NR) }'  

# Completion message
echo "Done."
date
echo ""