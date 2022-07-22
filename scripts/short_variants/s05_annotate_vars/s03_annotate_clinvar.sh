#!/bin/bash
# s04_annotate_clinvar - Annotate variants using ClinVar
# Anisha Thind, 4June2022

# Intended use:
# ./s04_annotate_clinvar.sh out_dir clinvar threads &> s04_annotate_clinvar.log
#  $1 (out_dir): output directory
#  $2 (clinvar): ClinVar VCF path
#  $3 (threads): number of threads


# ClinVar annotation of input VCF

# stop at runtime errors
set -euo pipefail

# starting message
echo -e "Script: s04_annotate_clinvar.sh\n"
date
echo ""

# set files and folders
out_dir="${1}"
vcf=$( find "${out_dir}" -name *.sorted.vcf.gz ) 
clinvar="${2}"
clinvar="${clinvar%.vcf.gz}.updated.vcf.gz"
threads="${3}"
annot_vcf="${vcf%%.*}.clinvar.vcf.gz"

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
elif [ ! -e "${clinvar}" ]; then
   echo "Error: ClinVar VCF file not found."
   exit 1
fi

# check update VCF exists
if [ -z "${vcf}" ]; then
  echo "Updated VCF file with IDs not found."
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

echo "Input VCF: ${vcf}"
echo "ClinVar VCF: ${clinvar}"
echo "Output VCF: ${annot_vcf}"
echo ""



# annotate variants using bcftools
echo "Annotating variants with ClinVar..."
bcftools annotate "${vcf}" \
   -a "${clinvar}" \
   -c "INFO" \
   --threads "${threads}" \
   -Oz \
   -o "${annot_vcf}"

# index VCF output file
bcftools index "${annot_vcf}"


echo "Number of INFO fields in input VCF:"
bcftools view -h "${vcf}" | grep "^##INFO" | wc -l
echo ""

echo "Number of INFO files in annotated VCF:"
bcftools view -h "${annot_vcf}" | grep "^##INFO" | wc -l
echo ""

echo "List of INFO fields in annotated VCF:"
echo ""
bcftools view -h "${annot_vcf}" | grep "^##INFO"
echo ""

echo "Number of variants with ClinVar annotations: "
bcftools query "${annot_vcf}" \
   -i 'ALLELEID != "."' \
   -f '%ID\tALLELEID\t%CLNSIG\t%CLNDN\n' | wc -l
echo ""

echo "Total number of variants in VCF: "
zgrep -v "^#" "${annot_vcf}" | wc -l
echo ""

echo "--- Variant counts by clinical significance ---"
bcftools query "${annot_vcf}" \
   -i 'ALLELEID != "."' \
   -f '%CLNSIG\n' \
   | sort \
   | uniq -c \
   | sort -r \
   | awk '{ printf("%s\t%-5s\n", $1, $2) }'

echo ""

# completion message
echo "Done."
date
echo ""
