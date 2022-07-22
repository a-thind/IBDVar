#!/bin/bash
# s00_start_annotation.sh - starts variant annotation pipeline

# Anisha Thind, 26Jun2022
# Parameters:
#   $1: (out_dir) output folder
#   $2: (clinvar) ClinVar VCF path
#   $3: (threads) Number of threads
#   $4: (log_dir) Log folder path

# stop at runtime errors
set -euo pipefail

# set parameters
out_dir="${1}"
clinvar="${2}"
vep="${3}"
cadd="${4}"
threads="${5}"
pipeline_log="${6}/s00_start_annotation.log"

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ !^[0-9]+ ]]; then
  threads=4
fi

echo ""
# check output dir
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi

out_dir="${out_dir}/s05_annotate_vars"
mkdir -p "${out_dir}"

echo -e "============================== Variant Annotation =============================\n" \
   |& tee "${pipeline_log}"
scripts_dir="short_variants/s05_annotate_vars"
echo -e "Data folder: ${out_dir}\n" \
   |& tee -a "${pipeline_log}"

echo -e "------------------------------ ClinVar Annotation -----------------------------\n" \
   |& tee -a "${pipeline_log}"
# check clinvar
"${scripts_dir}"/s01_check_clinvar_vcf.sh "${out_dir}" "${clinvar}" \
    "${threads}" \
   |& tee -a "${pipeline_log}"

"${scripts_dir}"/s02_rename_clinvar_chrs.sh "${out_dir}" "${clinvar}" \
    "${threads}" \
   |& tee -a "${pipeline_log}"

# # annotate variants using clinvar
# "${scripts_dir}"/s03_annotate_clinvar.sh "${out_dir}" "${clinvar}" \
#     "${threads}" \
#     |& tee -a "${pipeline_log}"

# echo -e "-------------------------------- VEP Annotation -------------------------------\n" \
#    |& tee -a "${pipeline_log}"
# "${scripts_dir}"/s04_annotate_vep.sh "${out_dir}" "${vep}" \
#    "${cadd}" \
#    "${threads}" \
#    |& tee -a "${pipeline_log}"

# "${scripts_dir}"/s05_split_vep_fields.sh "${out_dir}" \
#    |& tee -a "${pipeline_log}"

# completion message
echo -e "\nVariants annotated."
date
echo ""
