#!/bin/bash
# s02_add_recomb_dists.sh - insert recombination map into bim file (add recombination distances)
# Anisha Thind, 8Jun2022

# Intended use:
# ./s02_add_recomb_dists.sh out_dir ibis genetic_map threads &> s02_add_recomb_dists.log
# Paramters:
#   out_dir: output folder
#   ibis: IBIS folder
#   genetic_map: recombination genetic map

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# start message
printf "Script:\ts02_add_recomb_dists\n"
date
echo ""

# files and folders
out_dir="${1}/s07_select_haploblocks/plink"
ibis=$2
genetic_map=$3
add_map_plink="${ibis}/add-map-plink.pl"
plink_dataset="${out_dir}/autosomal_snps"


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