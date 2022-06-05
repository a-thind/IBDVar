#!/bin/bash

# Anisha Thind, 4June2022

# Intended use:
# ./s03_rename_clinvar_chrs.sh $VCF &> s03_rename_clinvar_chrs.log

# stop at runtime errors
set -e

# start message
echo $0
date
echo ""

# files and folders
VCF=$1
chr_map="/home/share/scripts/s04_annotate_vars/s02_clinvar_chr_map.txt"
clinvar="/home/share/resources/clinvar/clinvar_20220507.vcf.gz"
updated_clinvar="${clinvar%%.*}.updated.vcf.gz"

# progress report
bcftools --version
date
echo ""

echo "Input VCF: ${1}"
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
bcftools view -h "${VCF}" | grep "^##contig" | head -n 25
echo "..."
echo ""

echo "--- Chromosome Map File ---"
cat "${chr_map}"
echo ""

echo "Updating ClinVar VCF..."
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
