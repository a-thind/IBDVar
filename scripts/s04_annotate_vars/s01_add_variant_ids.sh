#!/bin/bash
# s01_add_variant_ids.sh - Adds ID for each variant
# Anisha Thind, 29May2022

# Intended use:
# ./s01_add_variant_ids.sh out_dir [threads] &> s01_add_variant_ids.log
# out_dir: output directory
# threads: number of threads (optional)



# stop at runtime errors, if any variable value is unset or if non-zero in pipe
set -euo pipefail

# Starting message
echo "Script: s01_add_variant_ids.sh"
date
echo ""

# Progress report
bcftools --version
date
echo ""

# set files and folders
out_dir="${1}"
threads="${2}"
data_dir="${out_dir}/s03_split_MA_sites"
out_dir="${out_dir}/s04_annotate_vars"
mkdir -p "${out_dir}"
vcf=$( find "${data_dir}" -name *.split_MA.vcf.gz ) 
filepath=$( basename "${vcf}" .vcf.gz ) 
basename="${filepath%%.*}"
out_vcf="${out_dir}/${basename}.ID.vcf.gz"

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ ^[0-9]+ ]]; then
  threads=4
fi

if [ -z "${vcf}" ]; then
  echo "Error: Missing MA split VCF file."
  exit 1
fi

# input VCF
echo ""
echo "Input VCF: ${vcf}"
echo "Output VCF (with IDs): ${out_vcf}"
echo ""

# Add variant IDs
echo "Adding variant IDs..."
bcftools annotate "${vcf}" \
  -I '%CHROM\_%POS\_%REF\_%ALT' \
  --threads "${threads}" \
  -Oz -o "${out_vcf}"

# index VCF
echo "Indexing VCF file..."
bcftools index "${out_vcf}"
echo ""

# Completion message
echo "Done."
date
echo ""
