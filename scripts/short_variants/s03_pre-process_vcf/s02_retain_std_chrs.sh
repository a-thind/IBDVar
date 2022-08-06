#!/bin/bash
# s02_retain_std_chrs.sh - retain standard chromosomes (1-22,X and Y) in annotated VCF
# Anisha Thind, 4Jun2022

# Intended use:
# ./s02_retain_std_chrs.sh out_dir threads &> s02_retain_std_chrs.sh
# out_dir: output directory
# threads: number of threads

# stop at runtime errors
set -euo pipefail

# start message
printf "Filtering Annotated VCF to retain only standard chromosomes\n\n"
printf "Script: s02_retain_std_chrs.sh\n"
date
echo ""

# files and folders
out_dir="${1}"
threads="$2"
vcf=$( find "${out_dir}" -name *.split_MA.vcf.gz ) 
std_chr="${vcf%.vcf.gz}.std_chr.vcf.gz"
updated_vcf="${std_chr%.std_chr*}.reheaded.vcf.gz"
header="${out_dir}/header.txt"

# check output
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi

# check update VCF exists
if [ -z "${vcf}" ]; then
  echo "Updated ClinVar VCF file not found."
  exit 1
fi

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ !^[0-9]+ ]]; then
  threads=4
fi

# progress report
bcftools --version
date
echo ""

echo "Input VCF file: ${vcf}"
echo "Output VCF file (with updated header): ${updated_vcf}"
echo ""

echo "Selecting variants from standard chromosomes only..."
bcftools annotate "${vcf}" \
   -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
   --threads "${threads}" \
   -Oz -o "${std_chr}"

# Index VCF
bcftools index "${std_chr}"
echo ""

# variant counts in original
echo "-- Variant counts in annotated VCF file --"
bcftools +counts "${vcf}"
echo ""
echo "-- Variant counts in standard chromosomes of annotated VCF --"
bcftools +counts "${std_chr}"
echo ""

# summary of variants per contig
echo "-- Variant counts in standard chromosomes --"
printf "%s\t%s\t\t%s\n" Contig Length Records
echo "------------------------------------"
bcftools index "${std_chr}" -s
echo ""

# update header of VCF
echo "Updating VCF header..."
bcftools view -h "${std_chr}"| zgrep -E '(^##contig=<ID=chr([0-9]{1,2}|X|Y),)|^##[^contig].*|^#CHROM.*' > "${header}"
# reheader VCF 
bcftools reheader "${std_chr}" \
   -h "${header}" \
   --threads "${threads}" \
   -o "${updated_vcf}"
echo ""

# index VCF
echo -e "Indexing VCF file...\n"
bcftools index "${updated_vcf}"

# verify changes have been made
echo "Number of contig fields in annotated VCF file:"
bcftools view -h "${std_chr}" | grep "^##contig" | wc -l
echo ""
echo "Number of contig fields in standard chromosome VCF file:"
bcftools view -h "${updated_vcf}" | grep "^##contig" | wc -l
echo ""

# Completion message
echo "Done."
date
echo ""
