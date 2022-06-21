#!/bin/bash
# short_variants.sh - Script for starting short variants prioritisation pipeline
# Anisha Thind, 15June2022

# stop at runtime errors
set -e
# stop pipeline if non-zero status error
set -o pipefail

# starting message
echo ""
echo "Short Variants Pipeline"
echo ""
date
echo ""


# Check input
# generate usage message
usage()
{
   printf "Usage: ./short_variants.sh [-c short_variants.config] [-m md5sum_file.md5sum]  \n\n" 
   printf "Options:\n"
   printf "\t-c,\tconfiguration file\n"
   printf "\t-m,\tmd5sum file for input VCF file\n"
   printf "\t-h,\tdisplay this help message\n"
   exit 0;
}
# display usage message if no args given
[ $# -eq 0 ] && usage

echo "----------------------------------- Settings ----------------------------------"
echo ""
# check parameter arguments
while getopts "c:m:h" arg; do
   case "${arg}" in
      h)
         usage
         exit 0
      ;;
      c)
         if [ -e "${OPTARG}" ]; then
            filename=` basename ${OPTARG} `
            if [ "${filename##*.}" == "config" ]; then
                echo "Configuration file: ${OPTARG}"
                echo ""
                . "${OPTARG}"
                printf 'Input VCF file:\t%-5s\n\n' "${IN_VCF}"
                printf 'Output folder:\t%-5s\n\n' "${OUT_DIR}"
                printf 'Number of threads:\t%-5s\n\n' "${THREADS}"
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
            filename=` basename "${OPTARG}"`
            if [ "${filename##*.}" == "md5sum" ]; then
               echo "--------------------------------- MD5SUM Check --------------------------------"
               echo ""
               echo "MD5SUM file: ${OPTARG}"
               echo ""
               # create folder for qc
               DATA_DIR="${OUT_DIR}/s01_short_vcf_qc"
               mkdir -p "${DATA_DIR}"
               LOG="${DATA_DIR}/s01_check_vcf.log"
               # run md5sum check
               . s01_short_vcf_qc/s01_check_vcf.sh "${IN_VCF}" "${OPTARG}" "${DATA_DIR}" &> "${LOG}"
               cat "${LOG}"
               echo ""
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

# Check VCF input
if [ ! -e "${IN_VCF}" ]
then
   echo "File ${IN_VCF} not found."
fi

#--------------------------------------------------------------------------------------------
# QC
#--------------------------------------------------------------------------------------------
echo "--------------------------- Computing BCFtools stats --------------------------"
echo ""
#LOG="${OUT_DIR}/s01_short_vcf_qc/s02_bcfstats.log"
#. s01_short_vcf_qc/s02_bcfstats.sh "${IN_VCF}" "${OUT_DIR}" &> "${LOG}"
#cat "${LOG}"
echo ""
echo ""
#----------------------------------------------------------------------------------------------
# Pre-processing
#----------------------------------------------------------------------------------------------
echo "--------------------------- Variant Pre-processing ----------------------------"
echo ""
echo "Filtering VCF..."
#DATA_DIR="${OUT_DIR}/s02_retain_pass_filter_vars"
#LOG="s01_retain_pass_filter_vars.log"
#. s02_retain_pass_filter_vars/s01_retain_pass_filter_vars.sh "${IN_VCF}" "${OUT_DIR}" &> "${LOG}"
#cat "${LOG}"
#mv "${LOG}" "${DATA_DIR}"
#echo ""
# plot stats
#LOG="s02_check_vcf_stats.log"
#. s02_retain_pass_filter_vars/s02_check_vcf_stats.sh "${OUT_DIR}" &> "${LOG}"
#cat "${LOG}"
#mv "${LOG}" "${DATA_DIR}"
echo ""
# multiallelic site parsing
echo "Parsing multiallelic sites..."
echo ""
DATA_DIR="${OUT_DIR}/s03_split_MA_sites"
LOG="s01_split_MA_sites.log"
. s03_split_MA_sites/s01_split_MA_sites.sh "${OUT_DIR}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""
echo "------------------------------ Variant Annotation -----------------------------"
echo ""

echo ""
# completion messages
echo "Done."
date
echo ""
