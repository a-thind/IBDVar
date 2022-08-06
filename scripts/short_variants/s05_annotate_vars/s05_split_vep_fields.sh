#!/bin/bash
# Anisha Thind, 10Jul2022
# s05_split_vep_fields.sh - splits VEP annotation fields into individual columns
# Parameters:
#   $1: (out_dir) output directory

# start message
echo "Script: s05_split_vep_fields.sh"
date
echo ""

# files and folders
out_dir="${1}"
in_vcf=$( find "${out_dir}" -name *filtered_ibd.sorted.clinvar.vep.vcf.gz )
tmp_vcf="${in_vcf%.vcf.gz}.tmp.vcf.gz"
out_vcf="${in_vcf%.vep.vcf.gz}.split-vep.vcf.gz"

# check output
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi

# check update VCF exists
if [ -z "${in_vcf}" ]; then
  echo "Updated ClinVar VCF file not found."
  exit 1
fi

# progress report
bcftools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo -e "Output split VCF file: ${out_vcf} \n"

echo -e "Parsing VEP annotation fields...\n"
bcftools +split-vep "${in_vcf}" \
    -p vep_ \
    -Oz \
    -o "${tmp_vcf}" \
    -c "-"


# index vcf
bcftools index "${tmp_vcf}"

echo -e "Removing unsplit CSQ field from VCF...\n"
bcftools annotate "${tmp_vcf}" \
    -x INFO/CSQ \
    --threads "${threads}" \
    -Oz \
    -o "${out_vcf}"

echo -e "Indexing split-VEP VCF...\n"
bcftools index "${out_vcf}"

# clean up temp files
rm "${tmp_vcf}"


echo ""
echo "Done."
date
echo ""
