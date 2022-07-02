#!/bin/bash
# s04_add_recomb_dists.sh - add recombination distances to filtered bim
# Anisha Thind, 8Jun2022

# Intended use:
# ./s04_add_recomb_dists.sh out_dir ibis genetic_map &> s04_add_recomb_dists.log
# Parameters:
#   $1: (out_dir) output folder
#   $2: (ibis) ibis folder path
#   $3: (genetic map) recombination genetic map for human genome

# stop script if non-zero exit status 
set -euo pipefail

echo -e "Script: s04_add_recomb_dists.sh\n"
date
echo ""

# files and folders
out_dir="${1}/s07_select_haploblocks/ibis"
ibis="${2}"
add_map_plink="${ibis}/add-map-plink.pl"
plink_dataset="${out_dir}/ibis"
genetic_map="${3}"
updated_plink="${out_dir}/ibis.recomb.bim"

if [ -z "${ibis}" ]; then
  echo "Error: Missing IBIS argument."
  exit 1
elif [ ! -d "${ibis}" ]; then
  echo "Error: IBIS path: ${ibis} is not a directory."
  exit 1  
fi

if [ ! -e "${genetic_map}" ]; then
  echo "Error: genetic map path does not exist."
  exit 1
fi

if [ ! -e "${add_map_plink}" ]; then
  echo "Error: cannot find add_map_plink.pl script."
  exit 1
fi

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