# s00_start_IBD_detection.sh - pipeline for selecting IBD regions.
# Anisha Thind, 26Jun2022

# Parameters:
#   $1: (out_dir) output directory
#   $2: (plink) plink path
#   $3: (threads) number of threads
#   $4: (ibis) IBIS folder path
#   $5: (genetic_map) genetic map path
#   $6: (ibis_mt) min number of markers for IBIS to call a segment 
#   $7: (truffle) truffle folder path
#   $8: (ibs1m) min number of markers for TRUFFLE to call a IBD1 type segment 
#   $9: (ibs2m) min number of markers for TRUFFLE to call a IBD2 type segment
#   $10: (phenogram) phenogram folder path
#   $11: (genome) Human genome text file with 3 columns (id, size, centromere) for phenogram
#   $12: pipeline log folder

# stop at runtime errors
set -euo pipefail

# set parameters
out_dir="${1}"
plink="${2}"
threads="${3}"
ibis="${4}"
genetic_map="${5}"
ibis_mt="${6}"
truffle="${7}"
ibs1m="${8}"
ibs2m="${9}"
genome="${10}"
log_dir="${11}"
pipeline_log="${log_dir}/s00_start_IBD_detection.log"

echo -e "================================ IBD DETECTION ================================\n" \
    |& tee "${pipeline_log}"

echo -e "Data folder: ${out_dir}\n" |& tee -a "${pipeline_log}"

echo -e "Script: s00_start_IBD_detection.sh\n" |& tee -a "${pipeline_log}"
date |& tee -a "${pipeline_log}"

# check output dir
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi

# check log dir
if [ -z "${log_dir}" ]; then
   echo "Error: Missing log directory argument."
   exit 1
elif [ ! -d "${log_dir}" ]; then
   echo "Error: log directory argument: ${log_dir} is not a directory."
   exit 1
fi

echo -e "\n--------------------------- Generating PLINK Dataset --------------------------\n"
scripts_dir="short_variants/s06_select_haploblocks"
"${scripts_dir}"/s01_make_plink_dataset.sh "${out_dir}" "${plink}" "${threads}" \
    |& tee -a "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n"

"${scripts_dir}"/s02_add_recomb_dists.sh "${out_dir}" "${ibis}" "${genetic_map}" \
    |& tee -a "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n"

"${scripts_dir}"/s03_retain_snps_with_recomb_dists.sh "${out_dir}" \
    "${plink}" \
    "${threads}" \
    |& tee -a "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n"

"${scripts_dir}"/s04_add_recomb_dists.sh "${out_dir}" \
    "${ibis}" \
    "${genetic_map}" \
    |& tee -a "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n"

"${scripts_dir}"/s05_select_haploblocks_ibis.sh "${out_dir}" \
    "${ibis}" \
    "${ibis_mt}" \
    "${threads}" \
    |& tee -a "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n"

"${scripts_dir}"/s06_select_haploblocks_truffle.sh "${out_dir}" \
    "${truffle}" \
    "${ibs1m}" \
    "${ibs2m}" \
    "${threads}" \
    |& tee -a "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n"

"${scripts_dir}"/s07_draw_ideogram.sh "${out_dir}" "${scripts_dir}" \
    |& tee -a "${pipeline_log}"

echo "IBD regions detected."
date
echo ""
echo -e "-------------------------------------------------------------------------------\n"
