#!/bin/bash
# s03_annotate_vep.sh - annotates variants using VEP
# Anisha Thind, 6Jul2022

# starting message
echo "Script: s03_annotate_vep.sh"
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data"
in_vcf="${data_dir}/sv/s03_annotate_sv"
vep="${base_dir}/tools/ensemble-vep/vep"
threads=4

# vep version
"${vep}" --version
date
echo ""

# run vep
"${vep}" -i "${in_vcf}" \
    --output_file "${out_vcf}" \
    --per_gene \
    --hgvs \
    --pubmed \
    --regulatory



# completions message
echo "Done."
date
echo ""