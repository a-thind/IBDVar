#!/bin/bash
# s03_retain_snps_with_recomb_dists.sh - retain only those SNPs that have recombination data
# Anisha Thind, 8Jun2022

# Intended use:
# ./s03_retain_snps_with_recomb_dists.sh out_dir plink threads &> s03_retain_snps_with_recomb_dists.log
# Parameters: 
#   out_dir: output folder
#   plink: Plink program folder
#   threads: number of threads

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# start message
printf "Script: s03_retain_snps_with_recomb_dists.sh\n\n"
date
echo ""

# files and folders
out_dir="${1}/s07_select_haploblocks"
plink="${2}/plink2"
threads=$3
plink_dataset="${out_dir}/plink/autosomal_snps"
# make dir for ibis data
out_dir="${out_dir}/ibis"
mkdir -p "${out_dir}"
filtered_plink="${out_dir}/ibis"
exclude_snps="${out_dir}/exclude_snps.txt"

# progress report
"${plink}" --version
date
echo ""
echo "Input plink dataset: ${plink_dataset}"
echo "Output filtered plink dataset: ${filtered_plink}"
echo "Output directory: ${out_dir}"
echo ""

echo "Excluding SNPs lacking recombination data..."
echo ""
# extract SNP IDs that lack recombination data (0 in third field)
awk '$3==0 {print $2}' "${plink_dataset}.recomb.bim" > "${exclude_snps}"

awk 'END{printf("Total number of SNPs: %s\n\n", NR)}' "${plink_dataset}.recomb.bim"

awk 'END{printf("Number of SNPS lacking recombination data: %s\n\n", NR)}' "${exclude_snps}"

"${plink}" -bfile "${plink_dataset}" \
    --exclude "${exclude_snps}" \
    --make-bed \
    --silent \
    --threads "${threads}" \
    --out "${filtered_plink}"
echo ""

# completion message
echo "Done."
date
echo ""
