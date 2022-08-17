#!/bin/bash
# s01_check_sv_vcf.sh - check
# Anisha Thind, 23Jun2022

# Example:
# ./s01_sv_vcf_stats.sh &> s01_sv_vcf_stats.log
# Parameters:
#   $1: (in_vcf) input SV VCF
#   $2: (out_dir) output folder 

# stop at runtime errors, any variable value is unset or non-zero status in pipe
set -euo pipefail

# starting message
echo "Script: s01_sv_vcf_stats.sh"
date
echo ""

# files and folders
in_vcf="${1}"
out_dir="${2}/s01_sv_vcf_qc"
basename=$( basename ${in_vcf} .sv.vcf.gz ) 

mkdir -p "${out_dir}"

# Progress report
bcftools --version
date
echo ""

echo ""
echo "Input VCF: ${in_vcf}"
echo "Output folder: ${out_dir}"
echo ""

echo "Indexing VCF..."
bcftools index -f "${in_vcf}"

echo -e "--- Summary of Variant Counts ---\n"
# extract SVTYPE info from VCF
bcftools query -f "%SVTYPE\n"  "${in_vcf}" \
    | sort \
    | uniq -c \
    | awk '{ printf("Number of %s:\t%s\n", $2, $1) }'

# number of SVs
bcftools view -H "${in_vcf}" | awk ' END{ printf("Number of SVs:\t%s\n\n", NR) } '

echo -e "--- Filters ---\n"
bcftools query -f "%FILTER\n" -i "FILTER='PASS'" "${in_vcf}" | awk 'END{ printf("Number of SVs passing all filters (PASS): %s\n\n", NR)}'

echo -e "Number of variants passing all filters by type and sample:\n"
bcftools query -f "[%SAMPLE %SVTYPE \n]" \
    -i "FILTER='PASS'" "${in_vcf}" \
    | sort \
    | uniq -c \
    | awk '{ printf("Number of %s (PASS) in %s: %s\n", $3, $2, $1) }'

echo -e "\n--- Variant counts by SV type and sample ---\n"
bcftools query -f "[%SAMPLE %SVTYPE \n]" -i "FILTER='PASS'" "${in_vcf}" \
    | sort \
    | uniq -c \
    | awk '{ printf("%s\t%s\t%s\n", $2, $3, $1) }' > "${out_dir}/sv_stats.tsv"

echo -e "Number of variants passing specific filters:\n"
bcftools query -f "%FILTER\n" "${in_vcf}" | sort | uniq -c | sort -nr 
echo ""

# Average length of insertions
bcftools query -f "%SVLEN\n" -i "SVTYPE='INS'" $in_vcf \
    | awk '{sum+=$1} END{printf("Average length of insertion: %s\n", sum / NR)}'
bcftools query -f "%SVLEN\n" -i "SVTYPE='DEL'" $in_vcf \
    | awk '{sum+=$1} END{printf("Average length for deletion: %s\n", sum / NR)}'
bcftools query -f "%SVLEN\n" -i "SVTYPE='DUP'" $in_vcf \
    | awk '{sum+=$1} END{printf("Average length for duplication: %s\n\n", sum / NR)}'

# Number of imprecise variants
echo "--- Imprecise Variants ---"
bcftools query -f "%SVTYPE %IMPRECISE\n" -i "IMPRECISE=1" "${in_vcf}" | sort | uniq -c | awk '{ printf("Number of imprecise %s: %s\n", $2, $1) }'
# total number of imprecise variants
bcftools query -f "%IMPRECISE\n" -i "IMPRECISE=1" "${in_vcf}" | awk 'END{ printf("Number of imprecise SVs: %s\n\n", NR)}'
echo ""

# Add Rscript
Rscript sv/s01_sv_vcf_qc/s02_sv_qc_stats.R "${in_vcf}" "${out_dir}"

# Completion message
echo "Done."
date
echo ""
