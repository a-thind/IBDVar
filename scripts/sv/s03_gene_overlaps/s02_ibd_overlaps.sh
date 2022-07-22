#!/bin/bash
# s02_ibd_overlaps.sh - compute overlaps wih IBD regions
# Anisha Thind, 6Jul2022

# Parameters:
#   $1: (out_dir) output directory

# start message
echo "Script: s02_ibd_overlaps.sh"
date
echo ""

# files and folders
out_dir="${1}"
vcf_dir="${out_dir}/s02_filter_sv"
in_vcf=$( find -name *.pass.read_support.vcf.gz )
ccds="${out_dir}/s03_gene_overlaps"

echo ""
echo "Input "
echo "Output"
echo ""

bedtools --version
date
echo ""

bcftools intersect 

echo "Done."
date
echo ""