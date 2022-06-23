#!/bin/bash
# s04_annotate_clinvar - Annotate variants using ClinVar
# Anisha Thind, 4June2022

# Intended use:
# ./s04_annotate_clinvar.sh out_dir clinvar threads &> s04_annotate_clinvar.log
# out_dir: output directory
# clinvar: ClinVar VCF path
# threads: number of threads (CPU)


# ClinVar annotation of input VCF

# stop at runtime errors
set -e
# stop if any variable value is unset
set -u
# stop pipeline if non-zero status
set -o pipefail

# starting message
echo "Variant Annotation using ClinVar"
printf "Script:\ts04_annotate_clinvar.sh\n"
date
echo ""

# set files and folders
out_dir=$1
vcf=` find "${out_dir}" -name *.ID.vcf.gz `
clinvar=$2
echo "clinvar"
clinvar="${clinvar%.vcf.gz}.updated.vcf.gz"
threads=$3
annot_vcf="${vcf%%.*}.clinvar.vcf.gz"

if [ -z "${vcf}" ]; then
  echo "VCF data file not found."
  exit 1
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
echo ""
bcftools view -h "${vcf}" | grep "^##INFO" | wc -l
echo ""

echo "Number of INFO files in annotated VCF:"
echo ""
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
   -f '%CLNSIG\n' | sort | uniq -c | sort -r | awk '{ printf("%s\t%-5s\n", $1, $2) }'
echo ""

# completion message
echo "Done."
date
echo ""
