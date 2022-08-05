#!/bin/bash
# s00_select_variants - selects impactful variants based on ClinVar and VEP annotations
# Anisha Thind, 25Jul2022

# Starting message
echo "Script: s00_select_variants"
date
echo ""

out_dir="${1}/s06_select_variants"
mkdir -p "${out_dir}"
MAF="${2}"
pipeline_log="${3}/s00_select_variants.log"
echo $out_dir

echo "=========== Selecting Variants ==============" > "${pipeline_log}"
short_variants/s06_select_variants/s01_filter_allele_freq.sh "${out_dir}" \
    "${MAF}" &>> "${pipeline_log}"

af_vcf=$( find "${out_dir}" -name *.AF.vcf.gz )

Rscript "${af_vcf}" "${out_dir}" &>> "${pipeline_log}"


echo "Done." >> "${pipeline_log}"
date >> "${pipeline_log}"
echo "" >> "${pipeline_log}"