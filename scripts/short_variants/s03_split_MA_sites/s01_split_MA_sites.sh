#!/bin/bash
# s01_split_MA_sites.sh - Splits alternate multiallele sites
# Anisha Thind, 18May2022

# Intended use:
# ./s01_split_MA_sites.sh in_dir out_dir &> s01_split_MA_sites.log
#   in_dir: input directory
#   out_dir: output directory


# stop at runtime errors
set -e
# stop if any variable value is unset
set -u 
# stop pipeline if non-zero status
set -o pipefail

# start message
echo -e "Multiallelic Site Parsing\n"
echo -e "Script: s01_split_MA_sites.sh\n"
date
echo ""

# set file and folder variables
in_dir="${1}"
out_dir="${2}"
vcf=$( find "${in_dir}" -name *.filtered.vcf.gz ) 

# Check source vcf was found
if [ -z "${vcf}" ]
then
    echo "Error: VCF file not found."
    exit 1
fi

# check output
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi

basename=$( basename "${vcf}" .vcf.gz ) 
ma_vcf="${out_dir}/${basename}.split_MA.vcf.gz"

# progress report
printf "Input VCF file: ${vcf}\n\n"
printf "Multi-allelic split VCF file: ${ma_vcf}\n\n"


# Make progress report
bcftools --version
date
echo ""

bcftools view -H "${vcf}" \
   | awk 'BEGIN{count=0} 
      $5~/,/{count++} 
      END{printf("Number of multi-allelic sites: %s\n\n", count)}'

# Split MA sites
printf "Split multiallelic sites...\n\n"
bcftools norm -m-any "${vcf}" -Oz -o "${ma_vcf}"
echo ""

# Indexing
printf "Indexing VCF...\n\n"
bcftools index "${ma_vcf}"
echo ""

# count variants
echo "Variant Counts"
echo "-----------------"
bcftools +counts "${ma_vcf}"
echo ""

# Completion message
echo "Done."
date
echo ""
