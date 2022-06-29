#!/bin/bash
# s00_start_pre-processing.sh - starts pass filtering pipeline and multi-allelic site parsing
# Parameters:
#  in_vcf: input vcf
#  out_dir: output folder

# Anisha Thind, 24Jun2022

# stop if runtime errors or non-zero status pipeline programs
set -euo pipefail

# parameters
in_vcf="${1}"
out_dir="${2}"

# check input vcf
if [ -z "${in_vcf}" ]; then
   echo "Error: Missing input VCF file argument." 
   exit 1
elif [ ! -e "${in_vcf}" ]; then
   echo "Error: VCF file: ${in_vcf} not found."
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

pipeline_log="${3}/s00_start_pre-processing.log"
data_dir="${out_dir}/s02_retain_pass_filter_vars"
mkdir -p "${data_dir}"
echo ""

echo -e "=========================== Variant Pre-processing ============================\n"
echo -e "---------------------- Filter Variants using Filter Field ---------------------\n"
base_dir="short_variants"
scripts_dir="${base_dir}/s02_retain_pass_filter_vars"
"${scripts_dir}"/s01_retain_pass_filter_vars.sh "${in_vcf}" \
   "${out_dir}" \
   |& tee -a "${pipeline_log}"

"${scripts_dir}"/s02_check_vcf_stats.sh "${out_dir}" \
    |& tee -a "${pipeline_log}"

# completion message
echo "Filtering using the Filter field completed."
date
echo ""

echo -e "------------------------- Multi-allelic Site Parsing --------------------------\n"
# set dirs
in_dir="${data_dir}"
out_dir="${out_dir}/s03_split_MA_sites"
mkdir -p "${out_dir}"
scripts_dir="${base_dir}/s03_split_MA_sites"

"${scripts_dir}"/s01_split_MA_sites.sh "${in_dir}" "${out_dir}" \
    |& tee -a "${pipeline_log}"

echo "Multi-allelic site parsing completed."
date
echo -e "\n-------------------------------------------------------------------------------\n"

echo "Variant pre-processing completed."
date
echo ""
