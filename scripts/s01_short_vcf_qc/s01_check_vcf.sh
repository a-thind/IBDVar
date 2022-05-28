#!/bin/bash

# s01_check_vcf.sh
# Anisha Thind, 10May2022

# Intended use:
# ./s01_check_vcf.sh vcf_file [md5_file] &> s01_check_vcf.log
# Performs md5sums check for a VCF file

# stop at runtime errors
set -e

# start message
echo "Check VCF file"
date
echo ""

# Store vcf filepath
VCF="${1}"
md5_file="${2}"

# if md5 file is provided
if [ ! -z "${md5_file}" ]
then
   if [ ! -e "${md5_file}" ]
   then
      cat "ERROR: File: ${md5_file} not found."
      cat "ERROR: computed checksums could not be verfied."
      exit 1
   else
      # check the md5sums match
      echo "Computed checksum:"
      md5sum "${VCF}" | grep -f <(awk '{print $1}' "${md5_file}")
      echo "md5sum file checksum:"
      cat "${md5_file}"
      echo ""

      if [ $? -eq 0 ]
      then
         echo "md5sum hashes match."
      else
         echo "WARNING: Computed checksums did NOT match."
      fi
   fi
   else
   # check m5sum of vcf file
   md5sum "${VCF}"
fi

echo ""

# end message
echo "Done."
date
