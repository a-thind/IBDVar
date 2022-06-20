#!/bin/bash
# s04_add_recomb_dists.sh - add recombination distances
# Anisha Thind, 8Jun2022

# Intended use:
# ./s04_add_recomb_dists.sh &> s04_add_recomb_dists.log

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s07_select_haploblocks/ibis"
ibis_dir="${base_dir}/tools/ibis"
add_map_plink="${ibis_dir}/add-map-plink.pl"
plink_dataset="${data_dir}/IHCAPX8_ibis.bim"
genetic_map="${base_dir}/resources/recomb-hg38/genetic_map_hg38_without_X.txt"
out_dataset="${data_dir}/IHCAPX8_ibis.recomb.bim"

# progress report
echo "Input plink dataset: ${plink_dataset}"
echo "Genetic map for hg38: ${genetic_map}"
echo "Output plink dataset: ${out_dataset}"
echo ""

echo "Inserting hg38 genetic map into bim file..."
echo ""
"${add_map_plink}" "${plink_dataset}" "${genetic_map}" > "${out_dataset}" 
echo ""

# completion message
echo "Done."
date
echo ""