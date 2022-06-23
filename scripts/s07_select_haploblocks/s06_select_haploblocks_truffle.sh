#!/bin/bash
# s06_select_haploblocks_truffle - identified haploblocks using truffle and annotated VCF as input

# Anisha Thind, 9Jun2022

# Intended use:
# ./s06_select_haploblocks_truffle.sh out_dir truffle ibs1m ibs2m threads &> s06_select_haploblocks_truffle.log
# Parameters:
#   out_dir: output folder
#   truffle: truffle folder path
#   threads: number of threads

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# starting message
printf "Script:\ts06_select_haploblocks_truffle\n"
date
echo ""

# files and folders
data_dir="${1}/s04_annotate_vars"
vcf=` find "${data_dir}" -name *.clinvar.reheaded.vcf.gz `
# create folder for truffle data
out_dir="${1}/s07_select_haploblocks/truffle"
mkdir -p "${out_dir}"
out_vcf="${out_dir}/truffle"
truffle="${2}/truffle"
ibs1m=$3
ibs2m=$4
threads=$5

if [ -z "${vcf}" ]; then
  echo "Annotated VCF file not found."
  exit 1
fi

# progress report
echo "Detecting haploblocks using truffle..."
echo ""
echo "Input VCF file (annotated with ClinVar): ${vcf}"
echo "Output folder: ${out_dir}"
echo ""

"${truffle}" --vcf "${vcf}" \
    --cpu "${threads}" \
    --segments \
    --ibs1markers "${ibs1m}" \
    --ibs2markers "${ibs2m}" \
    --out "${out_vcf}" 
echo ""

echo "Number of IBD segments identified:"
awk 'END{print NR-1}' "${out_vcf}.segments"
echo ""

# completion message
echo "Done."
date
echo ""