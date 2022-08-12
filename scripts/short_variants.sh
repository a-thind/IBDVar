#!/bin/bash
# short_variants.sh - Main script for starting short variants prioritisation pipeline
# Anisha Thind, 15June2022

# stop at runtime errors and if pipe contains a non-zero status
set -eo pipefail

# function to generate usage message
usage()
{
   echo -e "Usage: ./short_variants.sh -c short_variants.config [-m md5sum_file.md5sum]\n" 
   echo -e "Options:\n"
   echo -e "  -C, configuration file\n"
   echo -e "  -m, md5sum file for input VCF file\n"
   echo -e "  -h, display this help message\n"
   exit 0;
}
# display usage message if no args given
[ $# -eq 0 ] && usage

# check parameter arguments
while getopts "C:m:h" arg; do
   case "${arg}" in
      h)
         usage
         exit 0
      ;;
      C)
         if [ -e "${OPTARG}" ]; then
            filename=$( basename ${OPTARG} ) 
            if [ "${filename##*.}" == "config" ]; then
               config="${OPTARG}"
               # reads config file
               utils/read_sh_config.sh "${config}"
               . "${config}"
               # make log directory
               log_dir="${out_dir}/logs"
               mkdir -p "${log_dir}"
               pipeline_log="${log_dir}/pipeline.log"
            else
                echo "File ${OPTARG} is not a config file."
                exit 1
            fi
         else
            echo "File ${OPTARG} is not found."
            exit 1   
         fi
      ;;
      m)
         if [ -e "${OPTARG}" ]; then
            filename=$( basename "${OPTARG}" ) 
            if [ "${filename##*.}" == "md5sum" ]; then
               md5sum="${OPTARG}"
            else
               echo "File ${OPTARG} is not a md5sum file." 
               exit 1
            fi
         else
            echo "File ${OPTARG} is not found."
            exit 1
         fi
      ;;
      \?)
         echo "Invalid option: ${OPTARG}"
         usage
         exit 1
      ;;
   esac
done

# shift processed params so next param is $1
shift $((OPTIND-1))

# if there is a remaining param display usage
if [ ! -z "${1}" ]; then
   usage
fi

# starting message
echo -e "Short Variants Pipeline\n" &> "${pipeline_log}"
date &>> "${pipeline_log}"
echo -e "\n----------------------------------- Settings ----------------------------------\n" \
   | tee -a "${pipeline_log}"
echo -e "Configuration file: ${config}\n" \
   &>> "${pipeline_log}"
printf 'Input VCF file:\t%s\n\n' "${in_vcf}" \
   &>> "${pipeline_log}"
printf 'Output folder:\t%s\n\n' "${out_dir}" \
   &>> "${pipeline_log}"
printf 'Number of threads: %s\n\n' "${threads}" \
   &>> "${pipeline_log}"
echo -e "Log files are created in ${log_dir}\n\n" \
   &>> "${pipeline_log}"

echo -e "=================================== Pipeline ===================================\n" &>> "${pipeline_log}"
short_variants/s01_short_vcf_qc/s00_start_qc.sh "${in_vcf}" \
    "${out_dir}" \
    "${md5sum}" \
    "${log_dir}" \
    &>> "${pipeline_log}"

#--------------------------- Variant Pre-processing ----------------------------
short_variants/s02_filter_short_vars/s00_start_pre-processing.sh "${in_vcf}" \
    "${out_dir}" \
    "${GQ}" \
    "${DP}" \
    "${threads}" \
    "${log_dir}" \
    &>> "${pipeline_log}"

#-------------------------------- IBD Detection --------------------------------
short_variants/s04_select_ibd_variants/s00_start_IBD_detection.sh "${out_dir}" \
   "${plink}" \
   "${threads}" \
   "${ibis}" \
   "${genetic_map}" \
   "${ibis_mt1}" \
   "${ibis_mt2}" \
   "${truffle}" \
   "${ibs1m}" \
   "${ibs2m}" \
   "${genome}" \
   "${log_dir}" \
   &>> "${pipeline_log}" 

#------------------------------ Variant Annotation -----------------------------
short_variants/s05_annotate_vars/s00_start_annotation.sh "${out_dir}" \
    "${clinvar}" \
    "${vep}" \
    "${cadd}" \
    "${threads}" \
    "${log_dir}" \
    &>> "${pipeline_log}"

#------------------------------ Variant selection ------------------------------
short_variants/s06_select_variants/s00_select_variants.sh "${out_dir}" \
   "${MAF}" \
   "${log_dir}" \
   &>> "${pipeline_log}"

# completion messages
echo "Pipeline completed." &>> "${pipeline_log}"
date &>> "${pipeline_log}"
echo "" &>> "${pipeline_log}"
