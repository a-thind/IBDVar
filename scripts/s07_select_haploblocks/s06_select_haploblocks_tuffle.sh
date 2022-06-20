#!/bin/bash

# Anisha Thind, 9Jun2022

# Intended use:
# ./s07_select_haploblocks_truffle.sh &> s07_select_haploblocks_truffle.log

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# starting message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s04_annotate"
out_dir="${base_dir}/data/s07_select_haploblocks/truffle"
mkdir -p "${out_dir}"
out_vcf="${out_dir}/IHCAPX8_vcf_truffle"
out_bim="${out_dir}/IHCAPX8_bim_truffle"
vcf="${data_dir}/IHCAPX8_dragen_joint.clinvar.reheaded.vcf.gz"
plink_dataset="${out_dir%/truffle}/IHCAPX8_ibis"
truffle="${base_dir}/tools/truffle/truffle"
plink2="${base_dir}/tools/plink2/plink2"

# progress report
echo "Detecting haploblocks using truffle..."
echo ""
echo "Input VCF file (annotated with ClinVar): ${vcf}"
echo "Output folder: ${out_dir}"
echo ""

"${truffle}" --vcf "${vcf}" \
    --cpu 4 \
    --segments \
    --ibs1markers 4000 \
    --ibs2markers 500 \
    --out "${out_vcf}" 
echo ""

echo "Number of IBD segments identified:"
awk 'END{print NR-1}' "${out_vcf}.segments"
echo ""

# completion message
echo "Done."
date
echo ""