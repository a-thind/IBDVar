#!/bin/bash
# s01_check_vcf.sh - Checks md5sum matches between file and md5sum file
# Anisha Thind, 10May2022

# Intended use:
# ./s01_check_vcf.sh vcf_file in_vcf md5sum &> s01_check_vcf.log
# Options:
#  in_vcf: input VCF file
#  md5sum: md5sum file (ending with .md5sum extension)


# stop at runtime errors
set -e
# stop if variable value is unset
set -u
# stop pipeline if non-zero status
set -o pipefail

# start message
echo -e "Script: s01_check_vcf.sh\n"
date
echo ""

# Store vcf filepath
in_vcf="${1}"
md5_file="${2}"


echo "Input VCF file: ${in_vcf}"
echo "MD5SUM file: ${md5_file}"
echo ""

# check vcf file exists
if [ ! -e "${in_vcf}" ]; then
    echo "Input VCF file not found."
    exit 1
fi

# if md5 file is provided
if [ ! -z "${md5_file}" ]
then
   if [ ! -e "${md5_file}" ]
   then
      cat "Error: File: ${md5_file} not found."
      cat "Error: computed checksums could not be verfied."
      exit 1
   else
      # check the md5sums match
      echo "Computed checksum:"
      echo ""
      md5sum "${in_vcf}" | grep -f <(awk '{print $1}' "${md5_file}")
      echo ""
      echo "md5sum file checksum:"
      echo ""
      cat "${md5_file}"
      echo ""

      if [ $? -eq 0 ]
      then
         echo ""
         echo "md5sum hashes match."
      else
         echo "WARNING: Computed checksums did NOT match."
      fi
   fi
   else
   # check m5sum of vcf file
   md5sum "${in_vcf}"
fi

echo ""

# end message
echo "Done."
date
