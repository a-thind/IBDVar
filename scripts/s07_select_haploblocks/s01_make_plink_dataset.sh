#!/bin/bash
# s01_make_plink_dataset.sh - create plink dataset from annotated VCF
# Anisha Thind, 7Jun2022

# Intended use:
# ./s01_make_plink_dataset.sh out_dir plink threads &> s01_make_plink_dataset.log
# Parameters:
#   out_dir: output directory
#   plink: Plink path
#   threads: number of threads

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# start message
printf 'Script:\ts01_make_plink_dataset\n'
date
echo ""

# files and folders
out_dir=$1
plink="${2}/plink2"
threads=$3
vcf=` find "${out_dir}" -name *.clinvar.reheaded.vcf.gz `
out_dir="${out_dir}/s07_select_haploblocks/plink"
mkdir -p "${out_dir}"
plink_dataset="${out_dir}/autosomal_snps"


# progress report
"${plink}" --version
date
echo ""
echo "Input VCF: ${vcf}"
echo "Plink dataset: ${plink_dataset}"
echo ""

echo "Creating plink dataset..."
"${plink}" --vcf "${vcf}" \
    --vcf-half-call "missing" \
    --autosome \
    --snps-only \
    --make-bed \
    --silent \
    --threads "${threads}" \
    --out "${plink_dataset}"

echo ""

# Completion message
echo "Done."
date
echo ""