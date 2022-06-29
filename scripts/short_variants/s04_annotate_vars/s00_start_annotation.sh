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
threads="${3}"
pipeline_log="${4}/s00_start_annotation.log"

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

echo -e "============================== Variant Annotation =============================\n"
scripts_dir="short_variants/s04_annotate_vars"
# create variant ids
"${scripts_dir}"/s01_add_variant_ids.sh "${out_dir}" "${threads}" \
   |& tee -a "${pipeline_log}"

echo -e "------------------------------ ClinVar Annotation -----------------------------\n"
# check clinvar
"${scripts_dir}"/s02_check_clinvar_vcf.sh "${out_dir}" "${clinvar}" \
    "${threads}" \
   |& tee -a "${pipeline_log}"

"${scripts_dir}"/s03_rename_clinvar_chrs.sh "${out_dir}" "${clinvar}" \
    "${threads}" \
   |& tee "${pipeline_log}"

# annotate variants using clinvar
"${scripts_dir}"/s04_annotate_clinvar.sh "${out_dir}" "${clinvar}" \
    "${threads}" \
    |& tee "${pipeline_log}"

"${scripts_dir}"/s05_retain_std_chrs.sh "${out_dir}" "${clinvar}" \
   "${threads}" \
   |& tee "${pipeline_log}"

# completion message
echo -e "\nVariants annotated."
date
echo ""
