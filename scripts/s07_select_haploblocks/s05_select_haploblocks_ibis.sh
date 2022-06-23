#!/bin/bash
# s05_select_haploblocks_ibis.sh - runs IBIS to detect IBD segments
# Anisha Thind, 8Jun2022

# Intended use:
# ./s05_select_haploblocks_ibis.sh out_dir ibis ibis_mt threads &> s05_select_haploblocks_ibis.log
# Parameters:
#   out_dir: output folder
#   ibis: IBIS folder path
#   mt: min markers required for IBIS to call IBD 
#   threads: number of threads

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# starting message
printf "Script:\ts05_select_haploblocks_ibis.sh\n"
date
echo ""

# files and folders
out_dir="${1}/s07_select_haploblocks/ibis"
plink_dataset="${out_dir}/ibis"
ibd_segs="${plink_dataset}.seg"
ibis="${2}/ibis"
ibis_mt=$3
threads=$4

# progress report
echo "Input plink dataset: ${plink_dataset}"
echo "Output IBD segments file: ${ibd_segs}"
echo ""
# run IBIS
echo "Detecting IBD segments using IBIS..."
"${ibis}" -bfile "${plink_dataset}" \
    -ibd2 \
    -mt "${ibis_mt}" \
    -hbd \
    -t "${threads}" \
    -noFamID \
    -f "${plink_dataset}"  

# Insert header into segment file
sed -i '1i sample1  sample2  chrom    phys_start_pos phys_end_pos    IBD_type   genetic_start_pos genetic_end_pos genetic_seg_length  marker_count error_count error_density' "${ibd_segs}"

# completion message
echo "Done."
date
echo ""