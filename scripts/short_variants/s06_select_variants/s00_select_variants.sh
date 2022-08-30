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

echo "=========== Selecting Variants ==============" |& tee "${pipeline_log}"
short_variants/s06_select_variants/s01_filter_allele_freq.sh "${out_dir}" \
    "${MAF}" |& tee -a "${pipeline_log}"
# Use AF 
af_vcf=$( find "${out_dir}" -name *.AF.vcf.gz )

Rscript short_variants/s06_select_variants/s02_read_vcf_into_R.R "${af_vcf}" \
    "${out_dir}"  \
    |& tee -a "${pipeline_log}"


echo "Done." |& tee -a "${pipeline_log}"
date |& tee -a "${pipeline_log}"
echo "" |& tee -a "${pipeline_log}"