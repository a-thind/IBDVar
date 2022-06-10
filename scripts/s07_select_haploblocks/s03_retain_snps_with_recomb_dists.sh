#!/bin/bash

# Anisha Thind, 8Jun2022

# Intended use:
# ./s03_retain_snps_with_recomb_dists.sh &> s03_retain_snps_with_recomb_dists.log

# start message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s07_select_haploblocks"
plink_dataset="${data_dir}/IHCAPX8_autosomal_snps"
plink2="${base_dir}/tools/plink2/plink2"
output_dataset="${data_dir}/IHCAPX8_ibis"
exclude_snps="${data_dir}/exclude_snps.txt"

# progress report
"${plink2}" --version
date
echo ""
echo "Plink dataset: ${plink_dataset}"
echo "Filtered plink dataset: ${output_dataset}"
echo ""

echo "Excluding SNPs lacking recombination data..."
echo ""
# extract SNP IDs that lack recombination data (0 in third field)
awk '$3==0 {print $2}' "${plink_dataset}.recomb.bim" > "${exclude_snps}"
echo "Number of SNPs lacking recombination data:"
wc -l "${exclude_snps}" | awk '{print $1}'
echo ""

"${plink2}" -bfile "${plink_dataset}" \
    --exclude "${exclude_snps}" \
    --make-bed \
    --out "${output_dataset}"

echo ""
# completion message
echo "Done."
date
echo ""
