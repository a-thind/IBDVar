#!/bin/bash
# s02_check_clinvar_vcf.sh - checks chromosome notation in ClinVar and VCF file with IDs
# Anisha Thind, 4Jun2022

# Intended use:
# ./s02_check_clinvar_vcf.sh out_dir clinvar &> s02_check_clinvar_vcf.log
# out_dir: output folder
# clinvar: clinvar VCF path

# Stop at runtime errors
set -euo pipefail

# Starting message
echo "Script: s01_check_clinvar_vcf.sh"
date
echo ""

# files and folders
out_dir="${1}"
clinvar="${2}"
data_dir="${out_dir%%/s05_annotate_vars}/s04_select_ibd_variants"
vcf=$( find "${data_dir}" -name *.sorted.vcf.gz ) 
clinvar_chr="${out_dir}/clinvar_chr.txt"
vcf_chr="${out_dir}/vcf_chr.txt"
chr_map="${out_dir}/clinvar_chr_map.txt"

if [ -z "${vcf}" ]; then
  echo "Error: IBD filtered VCF file not found."
  exit 1
fi

# check output
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
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

# contig names from clinvar vcf
bcftools --version
date
echo ""

# generate progress report
echo "Input VCF: ${vcf}"
echo "ClinVar VCF: ${clinvar}"
echo ""

# Index clinvar vcf
echo "Indexing ClinVar VCF..."
bcftools index -f "${clinvar}"

echo "Checking ClinVar VCF..."
# for each vcf get contig names
bcftools index "${clinvar}" -s | awk 'NR<25{ print $1 }' > "${clinvar_chr}"
bcftools index "${vcf}" -s | awk 'NR<25{ print $1 }' > "${vcf_chr}"

echo "Creating translation table for chromosomes..."
# combine chromosome names for translation table text file
paste -d"\t" "${clinvar_chr}" "${vcf_chr}" > "${chr_map}"
echo "Chromosome translation table created."
# clean up contig names files
rm "${clinvar_chr}" "${vcf_chr}"
echo ""

echo "Reference in VCF data:"
echo ""
bcftools view -h "${vcf}" | grep "^##reference"
echo ""

echo "Reference in ClinVar:"
echo ""
bcftools view -h "${clinvar}" | grep "^##reference"
echo ""

# Completion message
echo "Done."
date
echo ""

