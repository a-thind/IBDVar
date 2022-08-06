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
# out_dir="${1}"
# vcf_dir="${out_dir}/s02_filter_sv"
# in_vcf=$( find -name *.pass.read_support.vcf.gz )
# ccds="${out_dir}/s03_gene_overlaps"
all_overlaps="/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/all_overlaps.bed"
ibd_bed="/home/share/data/output/IHCAPX8/s04_select_ibd_variants/ibis/ibd.bed"
ibd_genes="/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/ibd_filtered.bed"
in_vcf="$/home/share/data/outputIHCAPX8/sv/s02_filter_sv/IHCAPX8_SV_dragen_joint.sv.pass.read_support.vcf.gz"
ccds_filtered="/home/share/resources/ccds/ccds_filtered.bed"
ibd_overlaps="/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/ibd_overlaps.bed"
echo ""
echo "Input "
echo "Output"
echo ""

bedtools --version
date
echo ""

bedtools intersect -a "${ibd_bed}" \
    -b "${ccds_filtered}" \
    -wa \
    -wb > "${ibd_genes}"

echo "Number of genes within IBD regions:"
wc -l "${ibd_genes}"

echo "Done."
date
echo ""