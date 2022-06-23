#!/bin/bash
# s01_split_MA_sites.sh - Splits alternate multiallele sites
# s01_split_MA_sites.sh
# Anisha Thind, 18May2022

# Intended use:
# ./s01_split_MA_sites.sh out_dir &> s01_split_MA_sites.log
# out_dir: output directory


# stop at runtime errors
set -e
# stop if any variable value is unset
set -u 
# stop pipeline if non-zero status
set -o pipefail

# start message
printf "Multiallelic Site Parsing\n"
printf "./s01_split_MA_sites.sh\n"
date
echo ""

# set file and folder variables
out_dir=$1
data_dir="${out_dir}/s02_retain_pass_filter_vars"
vcf=` find "${data_dir}" -name *.pass_filtered.vcf.gz `
basename=` basename ${vcf} .vcf.gz `
out_dir="${out_dir}/s03_split_MA_sites"
mkdir -p "${out_dir}"
MA_VCF="${out_dir}/${basename}.split_MA.vcf.gz"

# progress report
printf "Input VCF file: ${vcf}\n\n"
printf "Multi-allelic split VCF file: ${MA_VCF}\n\n"
echo "${vcf}"
# Check source vcf was found
if [ -z "${vcf}" ]
then
    echo "Error: VCF file not found."
    exit 1
fi

# Make progress report
bcftools --version
date
echo ""

# Split MA sites
printf "Split multiallelic sites...\n\n"
bcftools norm -m-any "${vcf}" -Oz -o "${MA_VCF}"
echo ""

# Indexing
printf "Indexing VCF...\n\n"
bcftools index "${MA_VCF}"
echo ""

# count variants
echo "Variant Counts"
echo "-----------------"
bcftools +counts "${MA_VCF}"
echo ""

# Completion message
echo "Done."
date
echo ""
