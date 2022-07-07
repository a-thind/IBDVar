#!/bin/bash
# s01_annotate_sv - annotate SVs using ClinVar
# Anisha Thind, 4Jul2022

# stop at runtime errors
set -eou pipefail

# starting message
echo "Script: s01_annotate_sv.sh"
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/sv/s03_annotate_sv"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.pass.read_support_updated.vcf.gz"
clinvar="${base_dir}/resources/clinvar/sv/nstd102.GRCh38.variant_call.vcf.gz"
out_vcf="${in_vcf%.vcf*}.clinvar.vcf.gz"
threads=4

# annotate variants using ClinVar
bcftools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo "Output annotated VCF file: ${out_vcf}"
echo -e "ClinVar VCF: ${clinvar}\n\n"

echo -e "Annotating Variants with ClinVar...\n"
bcftools annotate "${in_vcf}" \
    -a "${clinvar}" \
    -c INFO \
    --threads "${threads}" \
    -Oz \
    -o "${out_vcf}"

echo -e "Indexing annotated VCF file...\n"
bcftools index "${out_vcf}"

echo "The number of variants lacking annotation:"
bcftools query -f "%INFO\n" -i 'DBVARID="."' "${out_vcf}" | wc -l
echo "Number of variants with ClinVar annotation:"
bcftools query -f "%INFO\n" -i 'DBVARID!="."' "${out_vcf}" | wc -l

# completion message
echo "Done."
date
echo ""