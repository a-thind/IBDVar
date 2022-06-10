#!/bin/bash

# Anisha Thind, 7Jun2022

# Intended use:
# ./s01_make_plink_dataset.sh &> s01_make_plink_dataset.log

# start message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s04_annotate"
in_vcf="${data_dir}/IHCAPX8_dragen_joint.clinvar.reheaded.vcf.gz"
out_dir="${base_dir}/data/s07_select_haploblocks"
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
    --double-id \
    --autosome \
    --snps-only \
    --make-bed \
    --out "${plink_dataset}"

echo ""

# Completion message
echo "Done."
date
echo ""