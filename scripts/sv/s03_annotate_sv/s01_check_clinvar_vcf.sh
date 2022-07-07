#!/bin/bash
# s01_check_clinvar_vcf.sh - Check ClinVar VCF file
# AT, 5Jul2022

# stop at runtime errors
set -eou pipefail

# start message
echo "Script: s01_check_clinvar_vcf.sh"
date
echo ""

# files and 
base_dir="/home/share"
data_dir="${base_dir}/data/sv/s02_filtering_sv"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.pass.read_support.vcf.gz"
clinvar="${base_dir}/resources/clinvar/sv/nstd102.GRCh38.variant_call.vcf.gz"
out_dir="${data_dir%/s02_*}/s03_annotate_sv"
basename=$( basename "${in_vcf}" .vcf.gz )
out_vcf="${out_dir}/${basename}_updated.vcf.gz"

# progress report
bcftools --version
date
echo ""

echo -e "Input VCF: ${in_vcf}\n"
echo -e "ClinVar VCF (SV): ${clinvar}\n\n"

# indexing clinvar vcf
bcftools index -f "${clinvar}"

# function to check chromosomes in ClinVar and input VCF files
# $1: ClinVar VCF file
# $2: Input VCF file
checkvcf() {
    echo "Checking chromosomes in ClinVar VCF:"
    bcftools query -f "%CHROM\n" "${1}" | sort | uniq -c | awk '{print $2}'

    echo -e "\nChecking chromosomes in Input VCF file:"
    bcftools query -f "%CHROM\n" "${2}" | sort | uniq -c | awk '{print $2}'
    echo ""
}

checkvcf "${clinvar}" "${in_vcf}"

echo -e "Ensuring chromosomes notation matches between ClinVar and input VCF files...\n"

# remove 'chr' prefix
zcat "${in_vcf}" \
    | sed -e 's/^chr//g' \
        -e 's/^##contig=<ID=chr/##contig=<ID=/g' \
        | bcftools view -Oz -o "${out_vcf}"

checkvcf "${clinvar}" "${out_vcf}"

echo -e "Indexing VCF file...\n"
bcftools index -f "${out_vcf}"


echo "Done."
date
echo ""