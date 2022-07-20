#!/bin/bash
# s01_filter_variants
# Anisha Thind, 7Jul2022

# stop at runtime errors
set -euo pipefail

# start message
echo "Script: s01_filter_vep.sh"
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/output/IHCAPX8/s04_annotate_vars"
in_vcf="${data_dir}/IHCAPX8_dragen_joint.clinvar.reheaded.split-vep.vcf.gz"
out_dir="${data_dir%/s04*}/s05_filter_vep_vars"
mkdir -p "${out_dir}"
out_vcf="${out_dir}/IHCAPX8_dragen_joint.clinvar.reheaded.VEP.AF.vcf.gz"
MAF=0.01

# progress report
bcftools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo -e "Output filtered VCF file: ${out_vcf}\n\n"

# filter variants
# select variants AF < 0.01 (rare)
echo -e "Selecting rare variants...\n"
bcftools view "${in_vcf}" \
    -i "vep_MAX_AF<${MAF} || vep_MAX_AF="."" \
    -Oz \
    -o "${out_vcf}"

zgrep -v "^#" "${out_vcf}" \
    | awk 'END{printf("Number of rare variants: %s\n\n", NR)}'



echo "Done."
date
echo ""