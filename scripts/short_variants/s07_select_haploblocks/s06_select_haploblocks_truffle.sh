#!/bin/bash
# s06_select_haploblocks_truffle - identified haploblocks using truffle and annotated VCF as input

# Anisha Thind, 9Jun2022

# Intended use:
# ./s06_select_haploblocks_truffle.sh out_dir truffle ibs1m ibs2m threads &> s06_select_haploblocks_truffle.log
# Parameters:
#   $1: (out_dir) output folder
#   $2: (truffle) truffle folder path
#   $3: (ibs1m) min number of markers for IBD1 type segment
#   $4: (ibs2m) min number of markers for IBD2 type segment
#   $5: (threads) number of threads

# stop script at runtime errors 
set -euo pipefail

# starting message
echo -e "Script: s06_select_haploblocks_truffle.sh\n"
date
echo ""

# files and folders
data_dir="${1}/s04_annotate_vars"
vcf=$( find "${data_dir}" -name *.clinvar.reheaded.vcf.gz ) 
# create folder for truffle data
out_dir="${1}/s07_select_haploblocks/truffle"
mkdir -p "${out_dir}"
out_vcf="${out_dir}/truffle"
truffle="${2}/truffle"
ibs1m="${3}"
ibs2m="${4}"
threads="${5}"

if [ -z "${vcf}" ]; then
  echo "Error: Annotated VCF file: ${vcf} not found."
  exit 1
fi

if [ ! -e "${truffle}" ]; then
  echo "Error: truffle path: ${truffle} does not exist."
  exit 1
fi

if [[ "${ibs1m}" =~ !^[0-9]+ ]]; then
  echo "Error: 'ibs1m' argument provided is non-numerical."
  exit 1
fi

if [[ "${ibs2m}" =~ !^[0-9]+ ]]; then
  echo "Error: 'ibs2m' argument provided is non-numerical."
  exit 1
fi

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ !^[0-9]+ ]]; then
  threads=4
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
awk 'END{ print NR-1 }' "${out_vcf}.segments"
echo ""

# completion message
echo "Done."
date
echo ""