#!/bin/bash

# Anisha Thind, 8Jun2022

# Intended use:
# ./s05_select_haploblocks_with_ibis.sh &> s05_select_haploblocks_with_ibis.log

# starting message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s07_select_haploblocks"
plink_dataset="${data_dir}/IHCAPX8_ibis"
ibis="${base_dir}/tools/ibis/ibis"

# progress report
echo "Plink dataset: ${plink_dataset}"
# run IBIS
"${ibis}" "${plink_dataset}.bed" \
    "${plink_dataset}.bim" \
    "${plink_dataset}.fam" \
    -f "${plink_dataset}"  


# completion message
echo "Done."
date
echo ""