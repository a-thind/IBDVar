#!/bin/bash
# s02_filter_read_support.sh - filter variants using read support evidence
# Anisha Thind, 4Jul2022

# starting message
echo "Script: s02_filter_read_support.sh"
date
echo ""

# files and folders
out_dir="${1}/s02_filter_sv"
in_vcf="${out_dir}/IHCAPX8_SV_dragen_joint.sv.pass.vcf.gz"
out_vcf="${in_vcf%.vcf.gz}.read_support.vcf.gz"
# filtering params
PR="${2}"
SR="${3}"


# progress report
bcftools --version
date
echo ""

echo ""
echo "Input VCF file: ${in_vcf}"
echo -e "Output VCF file: ${out_vcf}\n"
echo ""

echo -e "Filtering SVs with PR > ${PR} and SR > ${SR}...\n"

# imprecise variants don't have split read but are still included
bcftools filter -i "(PR>${PR} && SR>${SR}) || IMPRECISE=1" "${in_vcf}" -Oz  -o "${out_vcf}"
# count variants before and after filtering
zgrep -v "^#" "${in_vcf}" \
    | awk 'END{printf("Number of variants BEFORE filtering: %s\n", NR)}' 
zgrep -v "^#" "${out_vcf}" \
    |awk 'END{printf("Number of variants AFTER filtering: %s\n\n", NR)}'

echo "Indexing file..."
bcftools index "${out_vcf}"

# completion message
echo "Done."
date
echo ""