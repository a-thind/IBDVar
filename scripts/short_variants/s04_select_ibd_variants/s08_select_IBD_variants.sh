#!/bin/bash
# S08_select_IBD_variants.sh - select IBD variants from pre-process VCF
# Anisha Thind, 12Jun2022

# Example:
# ./s08_select_IBD_variants.sh &> s08_select_IBD_variants.log
# Parameters:
#   $1: (out_dir) output directory

# starting message
echo "Script: s08_select_IBD_variants.sh"
date
echo ""

# files and folders
out_dir="${1}"
ibd_dir="${out_dir}/ibis"
ibd_seg="${ibd_dir}/ibis.seg"
parent_dir=$( dirname "${out_dir}" )
vcf_dir="${parent_dir}/s03_pre-process_vcf"
in_vcf=$( find "${vcf_dir}" -name *.ID.vcf.gz )
ibd_bed="${ibd_dir}/ibd.bed"
basename=$( basename "${in_vcf}" .ID.vcf.gz )
out_vcf="${ibd_dir}/${basename}.filtered_ibd.vcf.gz"
sort_vcf="${out_vcf%.vcf.gz}.sorted.vcf.gz"

# progress report
bcftools --version
date
echo ""

# create bed file from segment file
awk -F"\t" '{OFS="\t"} NR>1{ print $3, $4, $5}' "${ibd_seg}" \
    | sed "s/^/chr/" > "${ibd_bed}"
sed 1i"CHROM\tPOS\tEND" "${ibd_bed}"

echo "Filtering variants in IBD regions..."
bcftools filter "${in_vcf}" \
    -R "${ibd_bed}" \
    --threads "${threads}" \
    -Oz \
    -o "${out_vcf}"

bcftools index -f "${out_vcf}"

echo "Sorting IBD VCF file..."
bcftools sort "${out_vcf}" \
    -Oz \
    -o "${sort_vcf}"

bcftools index -f "${sort_vcf}"

echo "Number of variants in IBD regions:"
zgrep -v "^#" "${sort_vcf}" | wc -l

# completion message
echo -e "\nDone."
date
echo ""