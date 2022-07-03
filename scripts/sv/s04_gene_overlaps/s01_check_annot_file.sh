#!/bin/bash
# s01_check_annot_file.sh - check annotation (CCDS) file chromosomes
# Anisha Thind, 25Jun2022

# stop at runtime errors, unset variable values and  
set -euo pipefail

echo -e "Script: s01_check_annot_file.sh\n"
date
echo ""

# set file and folders
base_dir="/home/share"
data_dir="${base_dir}/data/sv/s03_filter_imprecise_vars"
out_dir="${base_dir}/data/sv/s04_gene_overlaps"
mkdir -p "${out_dir}"
annot="${base_dir}/resources/ccds/CCDS.current.txt"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.pass.precise.vcf.gz"
out_vcf="${out_dir}/IHCAPX8_SV_dragen_joint.sv.pass.precise.vcf"

# progress report
bcftools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo "Annotation GTF: ${annot}"
echo ""
echo ""
echo "Ensuring names are consistent..."
echo "Chromosomes in VCF file: "
bcftools view -H "${in_vcf}"| cut -f1 | sort | uniq 
echo "Chromosomes in annotation file:"
cut -f1 "${annot}" | sort | uniq 


# prepare vcf file
gunzip -fk "${in_vcf}" 
sed "s/^chr//g" "${in_vcf%.gz}" >  "${out_vcf}"



# Completion message
echo "Done."
date
echo ""
