#!/bin/bash

# Anisha Thind, 6Jun2022

# Intended use:
# ./thin_vcf_test.sh $VCF &> thin_vcf_test.log

# stop at runtime errors
set -e

# Starting message
echo "Script: thin_vcf_test.sh"
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/output/IHCAPX8_1/s04_annotate_vars"
VCF="${data_dir}/IHCAPX8_dragen_joint.clinvar.std_chr.vcf.gz"
output_vcf="${data_dir}/thinned.recode.vcf"


# progress report
vcftools --version
date
echo ""

echo "Input VCF: ${VCF}"
echo "Output VCF: ${output_vcf}"
echo ""

# Thin VCF
echo "Thinning VCF..."
vcftools --gzvcf "${VCF}" \
   --thin 100000 \
   --recode \
   --recode-INFO-all \
   --out "${data_dir}/thinned"
echo ""

echo "Number of contigs in input VCF header:"
bcftools view -h "${VCF}" | grep "^##contig" | wc -l
echo ""
echo "Number of contigs in output VCF header:"
bcftools view -h "${output_vcf}" | grep "^##contig" | wc -l
echo ""

echo "Variant counts in input VCF:"
bcftools +counts "${VCF}"
echo ""
echo "Variant counts in output vcf:"
bcftools +counts "${output_vcf}"
echo ""

# Completion message
echo "Done."
date
echo ""