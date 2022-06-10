#!/bin/bash

# Anisha Thind, 4June2022

# Intended use:
# ./s04_annotate_clinvar.sh $VCF &> s04_annotate_clinvar.log

# ClinVar annotation of input VCF

# stop at runtime errors
set -e

# starting message
echo $0
date
echo ""

# set files and folders
VCF=$1
clinvar="/home/share/resources/clinvar/clinvar_20220507.updated.vcf.gz"
annot_VCF="${VCF%%.*}.clinvar.vcf.gz"

# progress report
bcftools --version
date
echo ""

echo "Input VCF: ${VCF}"
echo "ClinVar VCF: ${clinvar}"
echo ""

# annotate variants using bcftools
echo "Annotating variants with ClinVar..."
bcftools annotate "${VCF}" \
   -a "${clinvar}" \
   -c "INFO" \
   --threads 4 \
   -Oz \
   -o "${annot_VCF}"

# index VCF output file
bcftools index "${annot_VCF}"


echo "Number of INFO fields in input VCF:"
echo ""
bcftools view -h "${VCF}" | grep "^##INFO" | wc -l
echo ""

echo "Number of INFO files in annotated VCF:"
echo ""
bcftools view -h "${annot_VCF}" | grep "^##INFO" | wc -l
echo ""

echo "List of INFO fields in annotated VCF:"
echo ""
bcftools view -h "${annot_VCF}" | grep "^##INFO"
echo ""

echo "Number of variants with ClinVar AlleleIDs: "
bcftools query "${annot_VCF}" \
   -i 'ALLELEID != "."' \
   -f '%ID\tALLELEID\t%CLNSIG\t%CLNDN\n' | wc -l
echo ""

echo "--- Variant counts by clinical significance ---"
bcftools query "${annot_VCF}" \
   -i 'ALLELEID != "."' \
   -f '%CLNSIG\n' | sort | uniq -c | sort -r
echo ""

# completion message
echo "Done."
date
echo ""
