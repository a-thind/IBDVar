#!/bin/bash
# s01_filter_imprecise_vars.sh - Filter out imprecise SVs
# Anisha Thind, 25Jun2022

# stop if any runtime errors, unset variable values or non-zero exit status in pipe
set -euo pipefail

# files and folders
base_dir="/home/share/data/sv"
out_dir="${base_dir}/s03_filter_imprecise_vars"
mkdir -p "${out_dir}"
data_dir="${base_dir}/s02_retain_pass_filter_vars"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.pass.vcf.gz"
basename=$( basename "${in_vcf}" .vcf.gz )
out_vcf="${out_dir}/${basename}.precise.vcf.gz"
threads=4

echo -e "Script: s01_filter_imprecise_vars\n"
date
echo ""

# progress report
echo "Input VCF: ${in_vcf}"
echo "Output VCF: ${out_vcf}"
echo "Output folder: ${out_dir}"
echo ""

bcftools view -H "${in_vcf}" | awk 'END{ printf("Total number of variants BEFORE filtering: %s\n\n", NR) }'  

echo -e "Filtering out variants that are imprecise SVs..."
bcftools filter -Oz -i "IMPRECISE=0" "${in_vcf}" --threads "${threads}" > "${out_vcf}"
echo ""

echo -e "Indexing VCF file...\n"
bcftools index "${out_vcf}"

bcftools view -H "${out_vcf}" | awk 'END{ printf("Total number of variants AFTER filtering: %s\n\n", NR) }'  

# Completion message
echo "Done."
date
echo ""