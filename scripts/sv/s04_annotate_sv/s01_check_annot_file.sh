#!/bin/bash
# s01_check_annot_file.sh - check annotation (GTF/GFF) file chromosomes
# Anisha Thind, 25Jun2022

# stop at runtime errors, unset variable values and  
set -euo pipefail

echo -e "Script: s01_check_annot_file.sh\n"
date
echo ""

# set file and folders
base_dir="/home/share"
data_dir="${base_dir}/data/sv/s03_filter_imprecise_vars"
annot="${base_dir}/resources/refseq/hg38.ncbiRefSeq.gtf.gz"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.pass.precise.vcf.gz"

# progress report
bcftools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo "Annotation GTF: ${annot}"
echo ""

echo "Chromosomes in VCF file: "
bcftools index "${in_vcf}" -s | awk 'NR<25{ print $1 }'

#gunzip "${in_vcf}" | bcftools view -H | sed -i '/^chr/ !s/chr/'

# chromosome pattern: /(chr([0-9]{1,2}|[A-Z])|^[0-9]{1,2}|[A-Z])$/

# Completion message
echo "Done."
date
echo ""
