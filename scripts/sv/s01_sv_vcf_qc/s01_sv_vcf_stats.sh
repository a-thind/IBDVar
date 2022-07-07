#!/bin/bash
# s01_check_sv_vcf.sh - check
# Anisha Thind, 23Jun2022

# Intended use:
# ./s01_sv_vcf_stats.sh &> s01_sv_vcf_stats.log

# stop at runtime errors, any variable value is unset or non-zero status in pipe
set -euo pipefail

# starting message
printf "Script: s01_sv_vcf_stats.sh\n\n"
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s00_source_data/IHCAPX8/s02_structural_variants"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.vcf.gz"
csv="${data_dir}"
stats_dir="${base_dir}/data/sv/s01_sv_vcf_qc/bcfstats"
basename=$( basename ${in_vcf} .sv.vcf.gz ) 
stats_file="${stats_dir}/${basename}.vchk"


# Progress report
bcftools --version
date
echo ""

echo "Input VCF: ${in_vcf}"
echo "Output folder: ${stats_dir}"
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
    | awk '{ printf("%s\t%s\t%s\n", $2, $3, $1) }' > "${stats_dir}/sv_stats.tsv"

echo -e "Number of variants passing specific filters:\n"
bcftools query -f "%FILTER\n" "${in_vcf}" | sort | uniq -c | awk '{ printf("%-5s %s\n", $1, $2) }'
echo ""

# Number of structural variants
bcftools query -f "%SVTYPE %SVLEN\n" -i "SVTYPE='BND'" "${in_vcf}" | sort | uniq -c
bcftools query -f "%FILTER %SVTYPE %SVLEN\n" -i "SVTYPE='INS' && FILTER='PASS'" "${in_vcf}" | sort | uniq -c
# Average length of insertions
bcftools query -f "%SVLEN\n" -i "SVTYPE='INS'" $in_vcf | awk '{sum+=$1} END{print sum / NR}'
bcftools query -f "%SVLEN\n" -i "SVTYPE='INS'" $in_vcf | awk '{sum+=$1} END{print sum / NR}'

# Number of imprecise variants
echo "--- Imprecise Variants ---"
bcftools query -f "%SVTYPE %IMPRECISE\n" -i "IMPRECISE=1" "${in_vcf}" | sort | uniq -c | awk '{ printf("Number of imprecise %s: %s\n", $2, $1) }'
# total number of imprecise variants
bcftools query -f "%IMPRECISE\n" -i "IMPRECISE=1" "${in_vcf}" | awk 'END{ printf("Number of imprecise SVs: %s\n\n", NR)}'
echo ""

# Add Rscript

# Completion message
echo "Done."
date
echo ""
