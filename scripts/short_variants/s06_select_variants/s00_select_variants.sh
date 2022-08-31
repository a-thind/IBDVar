#!/bin/bash
# s00_select_variants - selects impactful variants based on ClinVar and VEP annotations
# Anisha Thind, 25Jul2022
# Options:
#    $1 - output directory (out_dir)
#    $2 - maximum allele frequency for rare variants in any of the populations from gnomAD, 1000 genomes project and ESP
#    $3 - log directory
#    $4 - (optional) list of genes of interest
# stop at runtime errors
set -eou pipefail

# Starting message
echo "Script: s00_select_variants"
date
echo ""

out_dir="${1}/s06_select_variants"
mkdir -p "${out_dir}"
MAF="${2}"
pipeline_log="${3}/s00_select_variants.log"
genes="${4}"
vars_out="${out_dir}/filtered_short_vars.tsv"

echo "=========== Selecting Variants ==============" |& tee "${pipeline_log}"
short_variants/s06_select_variants/s01_filter_allele_freq.sh "${out_dir}" \
    "${MAF}" |& tee -a "${pipeline_log}"
# Use AF 
af_vcf=$( find "${out_dir}" -name *.AF.vcf.gz )

Rscript short_variants/s06_select_variants/s02_read_vcf_into_R.R "${af_vcf}" \
    "${out_dir}"  \
    |& tee -a "${pipeline_log}"


if [ ! -z "${genes}" ]; then

    # check if file is excel spreadsheet
    if [[ "${genes}" == *.xlsx ]]; then
        in2csv "${genes}" > "${genes%.*}.csv"
        genes="${genes%.*}.csv"
    fi

    echo -e "Filtering by genes of interest...\n"
    grep -w -f "${genes}" "${vars_out}" > "${vars_out}"
fi

echo "Done." |& tee -a "${pipeline_log}"
date |& tee -a "${pipeline_log}"
echo "" |& tee -a "${pipeline_log}"