#!/bin/bash
# s01_check_vcf.sh Checks md5sum matches between file and md5sum file
# Anisha Thind, 10May2022

# Intended use:
# ./s01_check_vcf.sh vcf_file in_vcf md5sum out_dir &> s01_check_vcf.log
# in_vcf: input VCF file
# md5sum: md5sum file (ending with .md5sum extension)
# out_dir: output directory (for creating "s01_short_vcf_qc" directory)

# stop at runtime errors
set -e
# stop if variable value is unset
set -u
# stop pipeline if non-zero status
set -o pipefail

# start message
echo $0
date
echo ""

# Store vcf filepath
in_vcf=$1
md5_file=$2
out_dir=$3


echo "Input VCF file: ${in_vcf}"
echo "MD5SUM file: ${md5_file}"
echo ""

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
