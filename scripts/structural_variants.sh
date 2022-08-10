#!/bin/bash
# structural_variants.sh - pipeline for analysing structural variants
# Anisha Thind, 4July2022

# stop at runtime errors
set -eo pipefail

# function to generate usage message
usage()
{
   echo -e "Usage: ./structural_variants.sh -c short_variants.config [-m md5sum_file.md5sum]\n" 
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
               utils/read_sv_config.sh "${config}"
               . "${config}"
               # make log directory
               log_dir="${out_dir}/logs"
               mkdir -p "${log_dir}"
               pipeline_log="${log_dir}/sv_pipeline.log"
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


# create log file
pipeline_log="${out_dir}/logs/sv_pipeline.log"

echo "Structural Variants Pipeline" |& tee "${pipeline_log}"
date |& tee -a "${pipeline_log}"

echo -e "\n================================== Settings ===================================\n" \
    |& tee -a "${pipeline_log}"
echo "Input VCF: ${sv_vcf}" |& tee -a "${pipeline_log}"
echo "Output folder: ${out_dir}" |& tee -a "${pipeline_log}"
echo "Logs folder: ${out_dir}/logs" |& tee -a "${pipeline_log}"

echo -e "\n================================== Quality Control ===================================\n" \
    |& tee -a "${pipeline_log}"

sv/s01_sv_vcf_qc/s01_sv_vcf_stats.sh "${sv_vcf}" "${out_dir}" \
    |& tee -a "${pipeline_log}"

echo -e "\n============================== Filtering variants =============================\n" \
    |& tee -a "${pipeline_log}"

sv/s02_filter_sv/s01_retain_pass_filter_vars.sh "${sv_vcf}" "${out_dir}" \
    |& tee -a "${pipeline_log}"

sv/s02_filter_sv/s02_filter_read_support.sh "${out_dir}" "${PR}" "${SR}" \
    |& tee -a "${pipeline_log}"

echo -e "\n============================== Detect SV Overlaps =============================\n" \
    |& tee -a "${pipeline_log}"

sv/s03_gene_overlaps/s01_ibd_overlaps.sh "${out_dir}" \
    "${ccds}" \
    "${threads}" \
    "${ibd_seg}" \
    |& tee -a "${pipeline_log}"

# completion message
echo "Pipeline completed."
date
echo ""

