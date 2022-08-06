#!/bin/bash
# Anisha Thind, 6Jul2022

# start message
echo "Script: s01_check_cnv_vcf.sh"
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s00_source_data/IHCAPX8/s03_copy_number_variants"
out_dir="${base_dir}/data/cnv/s01_check_cnv_vcf"
mkdir -p "${out_dir}"
in_vcf="${data_dir}/IHCAPX8_CNV_dragen_joint.cnv.vcf.gz"
md5Sum="${in_vcf}.md5sum"

echo "Input VCF file: ${in_vcf}"

bcftools --version
date
echo ""

# check md5sums
../../utils/check_md5sum.sh "${in_vcf}" "${md5sum}"

bcftools query -f "%SVTYPE\n" "${in_vcf}" | sort | uniq -c

# Print CNV types stats
echo "--- Summary of Types of CNVs ---"
bcftools query -f "%ID\n" "${in_vcf}" | awk -F':' '{print $2}' | sort | uniq -c | sort -nr
# Completion message
echo "Done."
date
echo ""