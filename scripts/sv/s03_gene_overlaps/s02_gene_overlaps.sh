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
all_overlaps="${out_dir}/s03_gene_overlaps/all_overlaps.bed"

in_vcf=$( find  *.pass.read_support.vcf.gz)
ccds_filtered="/home/share/resources/ccds/ccds_filtered.bed"
ibd_overlaps="/home/share/data/output/IHCAPX8/sv/s03_gene_overlaps/ibd_overlaps.bed"
# gene list files
genes="${4}"
genes_base="${genes%%.*}"
genes_bed="${genes_base}.bed"
genes_bed_sort="${genes_base}_sort.bed"

echo ""
echo "Input "
echo "Output"
echo ""

bedtools --version
date
echo ""

echo "Input VCF: ${in_vcf}"
echo "Gene list file: ${genes}"
echo "Output Bed file: ${gene_overlaps}"

# check if file is excel spreadsheet
if [ "${genes}" == "*.xlsx" ]; then
    in2csv "${genes}" > "${genes##.}.csv"
    genes="${genes##.}.csv"
fi

# extract matching gene list lines from filtered_bed
echo -e "Finding overlaps...\n"
grep -w -f "${genes}" "${filtered_ccds}" > "${genes_bed}"

# sort genes bed file
sort -k1,1 -k2,2n -k3,3n --parallel "${threads}" "${genes_bed}" > "${genes_bed_sort}"

bedtools intersect -a "${ibd_bed}" \
    -b "${ccds_filtered}" \
    -wa \
    -wb > "${ibd_genes}"

# genes of interest
bedtools intersect -a "${in_vcf}" \
    -b "${genes_bed_sort}" \
    -wa \
    -wb > "${gene_overlaps}"

awk 'END{printf("Number of SV overlaps with genes of interest: %s\n\n", NR)}' "${gene_overlaps}"


echo "Number of genes within IBD regions:"
wc -l "${ibd_genes}"

echo "Done."
date
echo ""