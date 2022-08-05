#!/bin/bash
# structural_variants.sh - pipeline for analysing structural variants
# Anisha Thind, 4July2022

# stop at runtime errors
set -eo pipefail

# Options
sv_vcf="${1}"
out_dir="${2}"
PR=8
SR=0.15
ccds="/home/share/resources/ccds"
threads=4
genes="/home/share/resources/li_2022_gene_list.txt"
# create log file
pipeline_log="${out_dir}/logs/sv_pipeline.log"

# check parameters
# check output dir
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi


echo "Structural Variants Pipeline" |& tee "${pipeline_log}"
date |& tee -a "${pipeline_log}"

echo -e "\n================================== Settings ===================================\n" \
    |& tee -a "${pipeline_log}"
echo "Input VCF: ${sv_vcf}" |& tee -a "${pipeline_log}"
echo "Output folder: ${out_dir}" |& tee -a "${pipeline_log}"
echo "Logs folder: ${out_dir}/logs" |& tee -a "${pipeline_log}"

echo -e "\n================================== Quality Control ===================================\n" \
    |& tee -a "${pipeline_log}"

sv/s01_sv_vcf_qc/s01_sv_vcf_stats.sh "${sv_vcf}" "${out_dir}" \
    |& tee -a "${pipeline_log}"

echo -e "\n============================== Filtering variants =============================\n" \
    |& tee -a "${pipeline_log}"

sv/s02_filter_sv/s01_retain_pass_filter_vars.sh "${sv_vcf}" "${out_dir}" \
    |& tee -a "${pipeline_log}"

sv/s02_filter_sv/s02_filter_read_support.sh "${out_dir}" "${PR}" "${SR}" \
    |& tee -a "${pipeline_log}"

echo -e "\n============================== Detect SV Overlaps =============================\n" \
    |& tee -a "${pipeline_log}"

sv/s03_gene_overlaps/s01_gene_overlaps.sh "${out_dir}" \
    "${ccds}" \
    "${threads}" \
    "${genes}" \
    |& tee -a "${pipeline_log}"

# completion message
echo "Pipeline completed."
date
echo ""

