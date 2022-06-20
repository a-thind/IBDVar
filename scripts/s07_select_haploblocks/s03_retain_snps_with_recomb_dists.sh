#!/bin/bash
# s03_retain_snps_with_recomb_dists.sh - retain only those SNPs that have recombination data
# Anisha Thind, 8Jun2022

# Intended use:
# ./s03_retain_snps_with_recomb_dists.sh &> s03_retain_snps_with_recomb_dists.log

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
data_dir="${base_dir}/data/s07_select_haploblocks/plink"
out_dir="${data_dir%plink}/ibis"
mkdir -p "${out_dir}"
plink_dataset="${data_dir}/IHCAPX8_autosomal_snps"
plink2="${base_dir}/tools/plink2/plink2"
output_dataset="${out_dir}/IHCAPX8_ibis"
exclude_snps="${out_dir}/exclude_snps.txt"

# progress report
"${plink2}" --version
date
echo ""
echo "Input plink dataset: ${plink_dataset}"
echo "Output filtered plink dataset: ${output_dataset}"
echo "Output directory: ${data_dir}"
echo ""

echo "Excluding SNPs lacking recombination data..."
echo ""
# extract SNP IDs that lack recombination data (0 in third field)
awk '$3==0 {print $2}' "${plink_dataset}.recomb.bim" > "${exclude_snps}"
awk 'END{printf("Number of SNPS lacking recombination data: %s\n", NR)}' "${exclude_snps}"
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
