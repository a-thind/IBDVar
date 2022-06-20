#!/bin/bash
# s01_make_plink_dataset.sh - create plink dataset from annotated VCF
# Anisha Thind, 7Jun2022

# Intended use:
# ./s01_make_plink_dataset.sh &> s01_make_plink_dataset.log

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# start message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s04_annotate"
in_vcf="${data_dir}/IHCAPX8_dragen_joint.clinvar.reheaded.vcf.gz"
out_dir="${base_dir}/data/s07_select_haploblocks/plink"
mkdir -p "${out_dir}"
plink_dataset="${out_dir}/IHCAPX8_autosomal_snps"
plink2="${base_dir}/tools/plink2/plink2"

# progress report
"${plink2}" --version
date
echo ""
echo "Input VCF: ${in_vcf}"
echo "Plink dataset: ${plink_dataset}"
echo ""

echo "Creating plink dataset..."
"${plink2}" --vcf "${in_vcf}" \
    --vcf-half-call "missing" \
    --autosome \
    --snps-only \
    --make-bed \
    --out "${plink_dataset}"

echo ""

# Completion message
echo "Done."
date
echo ""