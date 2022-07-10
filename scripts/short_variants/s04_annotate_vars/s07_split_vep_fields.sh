#!/bin/bash
# Anisha Thind, 10Jul2022
# s07_split_vep_fields.sh - splits VEP annotation fields into individual columns

# start message
echo "Script: s07_split_vep_fields.sh"
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/output/IHCAPX8/s04_annotate_vars"
in_vcf="${data_dir}/IHCAPX8_dragen_joint.clinvar.reheaded.VEP.vcf.gz"
tmp_vcf="${in_vcf%.vcf.gz}.tmp.vcf.gz"
out_vcf="${in_vcf%.VEP*}.split-vep.vcf.gz"

# progress report
bcftools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo -e "Output split VCF file: ${out_vcf} \n"

echo -e "Parsing VEP annotation fields...\n"
bcftools +split-vep "${in_vcf}" \
    -c "-" \
    -p vep_ \
    -Oz \
    -o "${tmp_vcf}"

# index vcf
bcftools index "${tmp_vcf}"

echo -e "Removing unsplit CSQ field from VCF...\n"
bcftools annotate "${tmp_vcf}" \
    -x INFO/CSQ \
    --threads "${threads}" \
    -Oz \
    -o "${out_vcf}"

echo -e "Indexing split-VEP VCF...\n"
bcftools index "${out_vcf}"

# clean up temp files
rm "${tmp_vcf}*"


echo ""
echo "Done."
date
echo ""
