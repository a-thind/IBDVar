#!/bin/bash
# s01_add_variant_ids.sh - Adds ID for each variant
# Anisha Thind, 29May2022

# Intended use:
# ./s01_add_variant_ids.sh out_dir [threads] &> s01_add_variant_ids.log
# out_dir: output directory
# threads: number of threads (optional)



# stop at runtime errors
set -e
# stop if any variable value is unset
set -u
# stop pipeline if non-zero status
set -o pipefail

# Starting message
echo "Add Variant IDs"
date
echo ""

# Progress report
bcftools --version
date
echo ""

# set files and folders
out_dir=$1
threads=$2
data_dir="${out_dir}/s03_split_MA_sites"
out_dir="${out_dir}/s04_annotate_vars"
mkdir -p "${out_dir}"
vcf=` find "${data_dir}" -name *.split_MA.vcf.gz `
filepath=` basename "${vcf}" .vcf.gz `
basename="${filepath%%.*}"
out_vcf="${out_dir}/${basename}".ID.vcf.gz

# If no threads specified then default is 4
if [ -z "${threads}" ]; then
  threads=4
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
