#!/bin/bash
# s02_gene_overlaps.sh - filter variant list to include only those overlapping with genes of interest list
# Anisha Thind, 6Jul2022

# Parameters:
#   $1: (out_dir) output directory
#   $2: (genes) a list of genes of interest

# start message
echo "Script: s02_gene_overlaps.sh"
date
echo ""

# files and folders
out_dir="${1}/s03_gene_overlaps"
sv_ibd_filter=$( find "${out_dir}" -name ibd_annotated_sv.tsv )
gene_overlaps="${out_dir}/ibd_annotated_sv_genes.tsv"

if [ -z "${sv_ibd_filter}" ]; then
    echo "Error: Could not find filtered IBD SV tsv file."
    exit 1
fi

# gene list files
genes="${2}"

bedtools --version
date
echo ""

echo "Gene list file: ${genes}"
echo "Output variants overlapping gene list: ${gene_overlaps}"

# check if file is excel spreadsheet
if [[ "${genes}" == *.xlsx ]]; then
    in2csv "${genes}" > "${genes%.*}.csv"
    genes="${genes%.*}.csv"
fi

# extract matching gene list lines from filtered_bed
echo -e "Finding overlaps...\n"
awk 'NR==1{print}' "${sv_ibd_filter}" > "${gene_overlaps}"
grep -w -f "${genes}" "${sv_ibd_filter}" >> "${gene_overlaps}"

awk 'END{printf("Number of overlaps with genes of interest: %s\n\n", NR-1)}' "${gene_overlaps}"

echo "Done."
date
echo ""
