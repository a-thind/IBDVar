#!/bin/bash
# s01_filter_variants
# Anisha Thind, 7Jul2022

# stop at runtime errors
set -euo pipefail

# start message
echo "Script: s01_filter_allele_freq.sh"
date
echo ""

# files and folders
out_dir="${1}"
in_vcf=$( find "${out_dir}" -name *.sorted.split-vep.vcf.gz)
out_dir="${out_dir}/s06_select_variants"
mkdir -p "${out_dir}"
out_vcf="${out_dir}/IHCAPX8_dragen_joint.clinvar.reheaded.VEP.AF.vcf.gz"
MAF="${2}"

# progress report
bcftools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo -e "Output filtered VCF file: ${out_vcf}\n"

zgrep -v "^#" "${in_vcf}" \
    | awk 'END{printf("Number of VEP annotated variants: %s\n\n", NR)}'

# filter variants
# select variants AF < 0.01 (rare)
echo -e "Selecting rare variants (AF < ${MAF})...\n"
bcftools view "${in_vcf}" \
    -i "vep_MAX_AF< ${MAF} || vep_MAX_AF='.'" \
    -Oz \
    -o "${out_vcf}"

zgrep -v "^#" "${out_vcf}" \
    | awk 'END{printf("Number of rare variants: %s\n\n", NR)}'



echo "Done."
date
echo ""