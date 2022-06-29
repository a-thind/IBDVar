#!/bin/bash
# s03_rename_clinvar_chrs.sh - Renames ClinVar chromosomes so that they match input VCF file
# Anisha Thind, 4June2022

# Intended use:
# ./s03_rename_clinvar_chrs.sh out_dir clinvar &> s03_rename_clinvar_chrs.log
#   $1: (out_dir): output directory
#   $2: (clinvar): path for ClinVar VCF
#   $3: (threads) number of threads

# stop at runtime errors
set -euo pipefail

# start message
echo "Rename ClinVar Chromosomes"
date
echo ""

# files and folders
out_dir="${1}/s04_annotate_vars"
clinvar="${2}"
threads="${3}"
vcf=$( find "${out_dir}" -name *.ID.vcf.gz ) 
chr_map="${out_dir}/clinvar_chr_map.txt"
updated_clinvar="${clinvar%%.*}.updated.vcf.gz"

# progress report
bcftools --version
date
echo ""

if [ -z "${vcf}" ]; then
  echo "VCF file not found."
  exit 1
fi

# check clinvar VCF exists
if [ -z "${clinvar}" ]; then
   echo "Error: Missing ClinVar VCF file path."
   exit 1
elif [ ! -e "${clinvar}"]; then
   echo "Error: ClinVar VCF file not found."
   exit 1
fi

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ !^[0-9]+  ]]; then
  threads=4
fi

# check translation map exists
if [ ! -e "${chr_map}" ]; then
  echo "Error: ${chr_map} not found."
fi

echo "Input VCF: ${vcf}"
echo "ClinVar: ${clinvar}"
echo "Chromosome map file: ${chr_map}"
echo ""

echo "--- Contigs in ClinVar VCF header ---"
echo ""
bcftools view -h "${clinvar}" | grep "^##contig"
echo ""

echo "--- Number of records in ClinVar VCF ---"
echo ""
bcftools view -H "${clinvar}" | wc -l
echo ""

echo "--- Contigs in input VCF ---"
echo ""
bcftools view -h "${vcf}" | awk -F"\n" '$0 ~ /^##contig*/ { count++; if (count < 25) print $0}' 
printf "...\n\n"

echo "--- Chromosome Map File ---"
cat "${chr_map}"
echo ""

printf "Updating ClinVar VCF...\n\n"
# rename clinvar chromosomes
bcftools annotate "${clinvar}" \
   --rename-chrs "${chr_map}" \
   --threads "${threads}" \
   -Oz -o "${updated_clinvar}"

# index updated clinvar
bcftools index "${updated_clinvar}"

# Updated ClinVar checks
echo "--- Contigs in updated ClinVar VCF header ---"
echo ""
bcftools view -h "${updated_clinvar}" | grep "^##contig"
echo ""

echo "--- Number of records in updated ClinVar VCF ---"
bcftools view -H "${updated_clinvar}" | wc -l
echo ""

# Completion message
echo "Done."
date
echo ""
