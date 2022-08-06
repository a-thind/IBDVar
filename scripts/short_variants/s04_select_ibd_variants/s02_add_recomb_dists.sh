#!/bin/bash
# s02_add_recomb_dists.sh - insert recombination map into bim file (add recombination distances)
# Anisha Thind, 8Jun2022

# Intended use:
# ./s02_add_recomb_dists.sh out_dir ibis genetic_map threads &> s02_add_recomb_dists.log
# Paramters:
#   $1: (out_dir) output folder
#   $2: (ibis) IBIS folder
#   $3: (genetic_map) recombination genetic map

# stop at runtime errors 
set -euo pipefail

# start message
printf "Script: s02_add_recomb_dists.sh\n"
date
echo ""

# files and folders
out_dir="${1}/plink"
ibis="${2}"
genetic_map="${3}"
add_map_plink="${ibis}/add-map-plink.pl"
plink_dataset="${out_dir}/autosomal_snps_uniq_pos"

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

# Downloaded genetic map for Hg38 from:
# https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/
# using:
# wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
# removed X chromosome entries (coded as chr 23) using:
# awk "$1!=23" genetic_map_hg38_withX.txt > genetic_map_hg38_without_X.txt

# progress report
echo "Input plink dataset: ${plink_dataset}"
echo "Genetic map file: ${genetic_map}"
echo "Output bim: ${plink_dataset}.recomb.bim"
echo ""

echo "Inserting genetic map into bim file..."
"${add_map_plink}" "${plink_dataset}.bim" "${genetic_map}" > "${plink_dataset}.recomb.bim"

echo ""

# completion message
echo "Done."
date
echo ""