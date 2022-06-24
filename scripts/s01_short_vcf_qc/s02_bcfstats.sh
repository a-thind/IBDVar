#!/bin/bash
# s02_bcfstats.sh - Generate bcftool stats and QC plots
# Anisha Thind, 10May2022

# Intended use:
# ./s02_bcfstats.sh vcf_file out_dir &> s02_bcfstats.log
# Performs QC for VCF files

# stop at runtime errors
set -e
# stop if a variable value is unset
set -u
# stop pipeline if non-zero status
set -o pipefail

# start message
printf "Script:\ts02_bcfstats.sh\n"
date
echo ""

# set files and folder variables
in_vcf=$1
out_dir=$2
stats_dir="${out_dir}/s01_short_vcf_qc/bcfstats"

# check VCF file exists
if [ ! -e "${in_vcf}" ]
then
   echo "VCF file: ${in_vcf} not found."
   exit 1
fi

# create output files
base_name=` basename ${in_vcf} `
# make stats directory
mkdir -p "${stats_dir}"
stats_file="${stats_dir}/${base_name}.vchk"

# Generate progress report
bcftools --version
echo ""
echo "Input VCF file: ${in_vcf}"
echo "Output folder: ${stats_dir}"
echo ""

# reindex vcf file
echo "Indexing VCF file..."
bcftools index -f "${in_vcf}"
echo ""

# count variants
echo "Variant Counts"
echo "-----------------"
bcftools +counts "${in_vcf}"
echo ""

# generate BCFstats
echo "Calculating BCF stats..."
bcftools stats -s - "${in_vcf}" > "${stats_file}"
echo ""

# generate BCF stats plot
echo "Generating BCF stats plots..."
plot-vcfstats -s -p "${stats_dir}" "${stats_file}"
echo ""

# Completion message
echo "Done."
date
echo ""
