#!/bin/bash
# s02_ibd_overlaps.sh - compute overlaps wih IBD regions
# Anisha Thind, 6Jul2022

# start message
echo "Script: s02_ibd_overlaps.sh"
date
echo ""

# files and folders
in_vcf=$( find "${}" -name *.pass.read_support.vcf.gz )
ibd_seg="${}"


echo ""
echo "Input "
echo "Output"
echo ""

bedtools --version
date
echo ""

bcftools intersect \

echo "Done."
date
echo ""