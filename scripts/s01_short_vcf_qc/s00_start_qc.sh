#!/bin/bash
# start_qc.sh - starts qc pipeline for short variants
# Parameters:
#  in_vcf: input vcf
#  out_dir: output folder
#  md5sum: md5sum file (optional)

# Anisha Thind, 24Jun2022

# stop if runtime errors or non-zero status pipeline programs
set -euo pipefail

# set parameters
in_vcf="${1}"
out_dir="${2}"
md5sum="${3}"
pipeline_log="${4}/s00_start_qc.log"

# check input vcf
if [ -z "${in_vcf}" ]; then
   echo "Error: Missing input VCF file argument." 
   exit 1
elif [ ! -e "${in_vcf}" ]; then
   echo "Error: VCF file: ${in_vcf} not found."
   exit 1
fi

# check output dir
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi


data_dir="${out_dir}/s01_short_vcf_qc"
mkdir -p "${data_dir}"


echo -e "============================== Short Variants QC ==============================\n"

# md5sum check if file is present
if [ -e "${md5sum}" ]; then
   echo -e "--------------------------------- MD5SUM Check --------------------------------\n"
   s01_short_vcf_qc/s01_check_vcf.sh "${in_vcf}" "${md5sum}" |& tee -a "${pipeline_log}"
   echo ""
   elif [ ! -z "${md5sum}" ]; then
      echo "Error: file ${md5sum} not found."
      exit 1
fi

echo -e "--------------------------- Computing BCFtools stats --------------------------\n"
# create bcfstats
s01_short_vcf_qc/s02_bcfstats.sh "${in_vcf}" "${out_dir}" |& tee -a "${pipeline_log}"

echo "QC completed."
date
echo ""
