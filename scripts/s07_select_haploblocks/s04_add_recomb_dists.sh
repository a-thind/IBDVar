#!/bin/bash
# s04_add_recomb_dists.sh - add recombination distances to filtered bim
# Anisha Thind, 8Jun2022

# Intended use:
# ./s04_add_recomb_dists.sh out_dir ibis genetic_map &> s04_add_recomb_dists.log
# Parameters:
#   out_dir: output folder
#   ibis: ibis folder path
#   genetic map: recombination genetic map for human genome

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

printf "Script:\ts04_add_recomb_dists.sh\n"
date
echo ""

# files and folders
out_dir="${1}/s07_select_haploblocks/ibis"
ibis=$2
add_map_plink="${ibis}/add-map-plink.pl"
plink_dataset="${out_dir}/ibis"
genetic_map=$3
updated_plink="${out_dir}/ibis.recomb.bim"

# progress report
echo "Input plink dataset: ${plink_dataset}"
echo "Genetic map for hg38: ${genetic_map}"
echo "Output plink dataset: ${updated_plink}"
echo ""

echo "Inserting hg38 genetic map into bim file..."
echo ""
"${add_map_plink}" "${plink_dataset}.bim" "${genetic_map}" > "${updated_plink}" 
echo ""

awk 'END{printf("Total number of SNPs: %s\n\n", NR)}' "${plink_dataset}.recomb.bim"

# completion message
echo "Done."
date
echo ""