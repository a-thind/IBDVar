#!/bin/bash
# s05_select_haploblocks_ibis.sh - runs IBIS to detect IBD segments
# Anisha Thind, 8Jun2022

# Intended use:
# ./s05_select_haploblocks_ibis.sh out_dir ibis ibis_mt threads &> s05_select_haploblocks_ibis.log
# Parameters:
#   $1: (out_dir) output folder
#   $2: (ibis) IBIS folder path
#   $3: (ibis_mt1) min markers required for IBIS to call an IBD segment IBD1
#   $4: (ibis_mt2) min markers required for IBIS to call an IBD segment IBD2
#   $5: (threads) number of threads

# stop script if non-zero exit status 
set -euo pipefail

# starting message
echo -e "Script: s05_select_haploblocks_ibis.sh\n"
date
echo ""

# parameters
out_dir="${1}/ibis"
ibis="${2}/ibis"
ibis_mt1="${3}"
ibis_mt2="${4}"
threads="${5}"
# output files
plink_dataset="${out_dir}/ibis"
ibd_segs="${plink_dataset}.seg"

# check for ibis plink dataset
if [[ ! -e "${ibis}" ]]; then
  echo "Error: Plink dataset for IBIS not found."
  exit 1
fi

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ !^[0-9]+ ]]; then
  threads=4
fi

# check ibis min threshold parameter
if [[ "${ibis_mt1}" =~ !^[0-9]+ ]]; then
  echo "Error: 'ibis_mt1' argument provided is non-numerical."
  exit 1
fi

# check ibis min threshold parameter
if [[ "${ibis_mt2}" =~ !^[0-9]+ ]]; then
  echo "Error: 'ibis_mt2' argument provided is non-numerical."
  exit 1
fi

# progress report
echo "Input plink dataset: ${plink_dataset}"
echo "Output IBD segments file: ${ibd_segs}"
echo ""
# run IBIS
echo "Detecting IBD segments using IBIS..."
"${ibis}" -bfile "${plink_dataset}" \
    -ibd2 \
    -mt "${ibis_mt1}" \
    -mt2 "${ibis_mt2}" \
    -t "${threads}" \
    -noFamID \
    -f "${plink_dataset}"  

# Insert header into segment file
sed -i '1i sample1  sample2  chrom    phys_start_pos phys_end_pos    IBD_type   genetic_start_pos genetic_end_pos genetic_seg_length  marker_count error_count error_density' "${ibd_segs}"
echo ""

# completion message
echo "Done."
date
echo ""