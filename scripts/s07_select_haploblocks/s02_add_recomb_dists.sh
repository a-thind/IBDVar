#!/bin/bash

# Anisha Thind, 8Jun2022

# Intended use:
# ./s02_add_recomb_dists.sh &> s02_add_recomb_dists.log

# start message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s07_select_haploblocks"
add_map_plink="${base_dir}/tools/ibis/add-map-plink.pl"
plink_dataset="${data_dir}/IHCAPX8_autosomal_snps"
genetic_map="${base_dir}/resources/recomb-hg38/genetic_map_hg38_without_X.txt"

# Download genetic map for Hg38 from:
#https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/
# using:
# wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
# removed X chromosome entries (coded as chr 23) using:
# awk "$1!=23" genetic_map_hg38_withX.txt > genetic_map_hg38_without_X.txt

# progress report
echo "Plink dataset: ${plink_dataset}"
echo "IBIS folder: ${ibis_folder}"
echo "Genetic map file: ${genetic_map}"
echo ""

echo "Inserting genetic map into bim file..."
"${add_map_plink}" "${plink_dataset}.bim" "${genetic_map}" > "${plink_dataset}.recomb.bim"

echo ""

# completion message
echo "Done."
date
echo ""