#!/bin/bash
# s03_rename_clinvar_chrs.sh - Renames ClinVar chromosomes so that they match input VCF file
# Anisha Thind, 4June2022

# Intended use:
# ./s03_rename_clinvar_chrs.sh out_dir clinvar &> s03_rename_clinvar_chrs.log
# out_dir: output directory
# clinvar: path for ClinVar VCF

# stop at runtime errors
set -e
set -u 
set -o pipefail

# start message
echo "Rename ClinVar Chromosomes"
date
echo ""

# files and folders
out_dir=$1
out_dir="${out_dir}/s04_annotate_vars"
vcf=` find ${out_dir} -name *.ID.vcf.gz `
chr_map="${out_dir}/clinvar_chr_map.txt"
clinvar=$2
echo "clinvar"
updated_clinvar="${clinvar%%.*}.updated.vcf.gz"

# progress report
bcftools --version
date
echo ""

if [ -z "${vcf}" ]; then
  echo "VCF file not found."
  exit 1
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
   --threads 4 \
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
