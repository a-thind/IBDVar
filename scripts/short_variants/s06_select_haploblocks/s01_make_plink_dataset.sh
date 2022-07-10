#!/bin/bash
# s01_make_plink_dataset.sh - creates plink dataset from annotated VCF
# Anisha Thind, 7Jun2022

# Intended use:
# ./s01_make_plink_dataset.sh out_dir plink threads &> s01_make_plink_dataset.log
# Parameters:
#   $1: (out_dir) output directory
#   $2: (plink) Plink path
#   $3: (threads) number of threads

# stop at runtime errors
set -euo pipefail

# start message
echo -e 'Script: s01_make_plink_dataset\n'
date
echo ""

# files and folders
out_dir="${1}"
plink="${2}"

# check plink path exists
if [ -z "${plink}" ]; then
   echo "Error: Missing PLINK path."
   exit 1
elif [ ! -d "${plink}" ]; then
   echo "Error: plink path argument: ${plink} is not a directory."
   exit 1
fi

plink="${plink}/plink2"

if [ ! -e "${plink}" ]; then
   echo "Error: PLINK executable not found."
   exit 1
fi

threads="${3}"
vcf=$( find "${out_dir}" -name *.clinvar.reheaded.vcf.gz ) 
out_dir="${out_dir}/s07_select_haploblocks/plink"
mkdir -p "${out_dir}"
plink_dataset="${out_dir}/autosomal_snps"
pos_pattern="${out_dir}/pos_patterns.txt"
dup_pos_snps="${out_dir}/dup_pos_snps.txt"
plink_uniq_pos="${plink_dataset}_uniq_pos"

# check output dir
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi



if [ -z "${vcf}" ]; then
   echo "ClinVar annotated VCF file not found."
fi

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ ^[0-9]+ ]]; then
  threads=4
fi

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
   --mind 1 \
   --geno 0.01 \
   --maf 0.05 \
   --make-bed \
   --silent \
   --threads "${threads}" \
   --out "${plink_dataset}"

echo ""

# remove SNPs with identical physical positions
# find duplicate positions and create patterns for 
# checking
cut -f2 "${plink_dataset}.bim" | sed 's/....$/_._./' \
   | sort \
   | uniq -d > "${pos_pattern}"

# find the duplicate patterns in the bim file
grep -o -f "${pos_pattern}" "${plink_dataset}.bim" > "${dup_pos_snps}"

"${plink}" -bfile "${plink_dataset}" \
   --exclude "${dup_pos_snps}" \
   --make-bed \
   --threads "${threads}" \
   --out "${plink_uniq_pos}"
echo ""

# Completion message
echo "Done."
date
echo ""