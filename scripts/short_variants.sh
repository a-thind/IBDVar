#!/bin/bash
# short_variants.sh - Script for starting short variants prioritisation pipeline
# Anisha Thind, 15June2022

# stop at runtime errors and if pipe contains a non-zero status
set -eo pipefail

# starting message
echo -e "Short Variants Pipeline\n"
echo ""
date
echo ""

# Check input
# generate usage message
usage()
{
   echo -e "Usage: ./short_variants.sh -c short_variants.config [-m md5sum_file.md5sum]\n" 
   echo -e "Options:\n"
   echo -e "\t-C, configuration file\n"
   echo -e "\t-m, md5sum file for input VCF file\n"
   echo -e "\t-h, display this help message\n"
   exit 0;
}
# display usage message if no args given
[ $# -eq 0 ] && usage

echo -e "----------------------------------- Settings ----------------------------------\n"
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
                echo "Configuration file: ${OPTARG}"
                echo ""
                config="${OPTARG}"
                utils/read_config.sh "${config}"
                . "${config}"
                printf 'Input VCF file:\t%s\n\n' "${in_vcf}"
                printf 'Output folder:\t%s\n\n' "${out_dir}"
                printf 'Number of threads:\t%s\n\n' "${threads}"
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


echo -e "=================================== Pipeline ===================================\n"
# make log directory
log_dir="${out_dir}/logs"
mkdir -p "${log_dir}"
pipeline_log="${log_dir}/pipeline.log"

short_variants/s01_short_vcf_qc/s00_start_qc.sh "${in_vcf}" \
   "${out_dir}" \
   "${md5sum}" \
   "${log_dir}" \
   |& tee -a "${pipeline_log}"

#--------------------------- Variant Pre-processing ----------------------------
short_variants/s02_retain_pass_filter_vars/s00_start_pre-processing.sh "${in_vcf}" \
   "${out_dir}" \
   "${log_dir}" \
   |& tee -a "${pipeline_log}"

# ------------------------------ Variant Annotation -----------------------------
short_variants/s04_annotate_vars/s00_start_annotation.sh "${out_dir}" \
   "${clinvar}" \
   "${threads}" \
   "${log_dir}" \
   |& tee -a "${pipeline_log}"

#-------------------------------- IBD Detection --------------------------------
short_variants/s07_select_haploblocks/s00_start_IBD_detection.sh "${out_dir}" \
   "${plink}" \
   "${threads}" \
   "${ibis}" \
   "${genetic_map}" \
   "${ibis_mt}" \
   "${truffle}" \
   "${ibs1m}" \
   "${ibs2m}" \
   "${phenogram}" \
   "${genome}" \
   "${log_dir}" \
   |& tee -a "${pipeline_log}" 


# completion messages
echo "Pipeline completed."
date
echo ""