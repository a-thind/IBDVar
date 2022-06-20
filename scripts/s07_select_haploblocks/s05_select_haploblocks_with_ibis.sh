#!/bin/bash
# s05_select_haploblocks_with_ibis.sh - runs IBIS to detect IBD segments

# Anisha Thind, 8Jun2022

# Intended use:
# ./s05_select_haploblocks_with_ibis.sh &> s05_select_haploblocks_with_ibis.log

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# starting message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s07_select_haploblocks/ibis"
plink_dataset="${data_dir}/IHCAPX8_ibis"
ibd_segs="${plink_dataset}.seg"
ibis="${base_dir}/tools/ibis/ibis"

# progress report
echo "Input plink dataset: ${plink_dataset}"
echo "Output IBD segments file: ${ibd_segs}"
echo ""
# run IBIS
echo "Detecting IBD segments using IBIS..."
"${ibis}" -bfile "${plink_dataset}" \
    -ibd2 \
    -mt 100 \
    -hbd \
    -t 4 \
    -noFamID \
    -f "${plink_dataset}"  

# Insert header into segment file
sed -i '1i sample1  sample2  chrom    phys_start_pos phys_end_pos    IBD_type   genetic_start_pos genetic_end_pos genetic_seg_length  marker_count error_count error_density' "${ibd_segs}"

# completion message
echo "Done."
date
echo ""