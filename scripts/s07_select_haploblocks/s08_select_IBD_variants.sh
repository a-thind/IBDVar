#!/bin/bash

# Anisha Thind, 12Jun2022

# Intended use:
# ./s08_select_IBD_variants.sh &> s08_select_IBD_variants.log

# starting message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data"
truffle_dir="${data_dir}/s07_select_haploblocks/truffle"
vcf_dir="${data_dir}/s04_annotate"
in_file="${data_dir}/IHCAPX8_vcf_truffle.segments"
in_vcf="${vcf_dir}/"
out_vcf="${in_file%.segments}_IBD"

# progress report
bcftools --version
date
echo ""

# find segments on the same chromosome
awk 'a[$4]++{count++}END{for (i in a)if(a[i]>1)print a[$1];}' "${in_file}" 
samples=`bcftools query -l "${in_vcf}"`


# completion message
echo "Done."
date
echo ""