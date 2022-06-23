#!/bin/bash
# s02_check_clinvar_vcf.sh - checks chromosome notation in ClinVar and VCF file with IDs
# Anisha Thind, 4Jun2022

# Intended use:
# ./s02_check_clinvar_vcf.sh out_dir clinvar &> s02_check_clinvar_vcf.log
# out_dir: output folder
# clinvar: clinvar VCF path

# Stop at runtime errors
set -e
# stop if any variable value is unset
set -u
# stop pipeline if any program has non-zero status
set -o pipefail

# Starting message
echo "Check ClinVar Chromosome Notation"
echo "./s02_check_clinvar_vcf.sh"
date
echo ""

# files and folders
out_dir=$1
out_dir="${out_dir}/s04_annotate_vars"
clinvar=$2
echo "${clinvar}"
vcf=` find "${out_dir}" -name *.ID.vcf.gz `
clinvar_chr="${out_dir}/clinvar_chr.txt"
vcf_chr="${out_dir}/vcf_chr.txt"
chr_map="${out_dir}/clinvar_chr_map.txt"

if [ -z "${vcf}" ]; then
  echo "VCF file not found."
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

