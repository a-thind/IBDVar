# s00_start_IBD_detection.sh - pipeline for selecting IBD regions.
# Anisha Thind, 26Jun2022

# Parameters:
#   $1: (out_dir) output directory
#   $2: (plink) plink path
#   $3: (mind) 
#   $4: (geno)
#   $5: (MAF)
#   $6: (threads) number of threads
#   $7: (ibis) IBIS folder path
#   $8: (genetic_map) genetic map path
#   $9: pipeline log folder

# stop at runtime errors
set -eou pipefail

# set parameters
out_dir="${1}"
plink="${2}"
mind="${3}"
geno="${4}"
MAF="${5}"
threads="${6}"
ibis="${7}"
genetic_map="${8}"
ibis_mt1="${9}"
ibis_mt2="${10}"
log_dir="${11}"
pipeline_log="${log_dir}/s00_start_IBD_detection.log"

echo -e "================================ IBD DETECTION ================================\n" \
    &> "${pipeline_log}"

echo -e "Data folder: ${out_dir}\n" &>> "${pipeline_log}"

echo -e "Script: s00_start_IBD_detection.sh\n" &>> "${pipeline_log}"
date &>> "${pipeline_log}"

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

out_dir="${out_dir}/s04_select_ibd_variants"
mkdir -p "${out_dir}"

echo -e "\n--------------------------- Generating PLINK Dataset --------------------------\n" \
    &>> "${pipeline_log}"
scripts_dir="short_variants/s04_select_ibd_variants"

echo -e "Data folder: ${out_dir}\n" &>> "${pipeline_log}"
"${scripts_dir}"/s01_make_plink_dataset.sh "${out_dir}" "${plink}" \
    "${threads}" \
    "${mind}" \
    "${geno}" \
    "${MAF}" \
    &>> "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n" \
    &>> "${pipeline_log}"

"${scripts_dir}"/s02_add_recomb_dists.sh "${out_dir}" "${ibis}" "${genetic_map}" \
    &>> "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n" \
    &>> "${pipeline_log}"

"${scripts_dir}"/s03_retain_snps_with_recomb_dists.sh "${out_dir}" \
    "${plink}" \
    "${threads}" \
    &>> "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n" \
    &>> "${pipeline_log}"

"${scripts_dir}"/s04_add_recomb_dists.sh "${out_dir}" \
    "${ibis}" \
    "${genetic_map}" \
    &>> "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n" \
    &>> "${pipeline_log}"

"${scripts_dir}"/s05_select_haploblocks_ibis.sh "${out_dir}" \
    "${ibis}" \
    "${ibis_mt1}" \
    "${ibis_mt2}" \
    "${threads}" \
    &>> "${pipeline_log}"

echo -e "-------------------------------------------------------------------------------\n" \
    &>> "${pipeline_log}"

"${scripts_dir}"/s07_draw_ideogram.sh "${out_dir}" "${scripts_dir}" \
    &>> "${pipeline_log}"

"${scripts_dir}"/s08_select_IBD_variants.sh "${out_dir}" \
    &>> "${pipeline_log}"

echo "IBD regions detected." &>> "${pipeline_log}"
date &>> "${pipeline_log}"
echo "" &>> "${pipeline_log}"
