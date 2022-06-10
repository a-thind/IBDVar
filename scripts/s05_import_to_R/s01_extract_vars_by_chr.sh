#!/bin/bash

# Anisha Thind, 10June2022

# Intended use:
# ./s01_extract_vars_by_chr.sh &> s01_extract_vars_by_chr.sh

# starting message
echo $0
date 
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s04_annotate"
out_dir="${base_dir}/data/s05_import_to_R"
mkdir -p "${out_dir}"
in_vcf="${data_dir}/IHCAPX8_dragen_joint.clinvar.reheaded.vcf.gz"
basename=` basename ${in_vcf} .vcf.gz `
out_prefix="${out_dir}/${basename}"

# progress report
bcftools --version
date 
echo ""

echo "Input VCF: ${in_vcf}"
echo "Output folder: ${out_dir}"
echo ""

echo "Extracting variants by chromosomes into individual VCF files..."
echo ""
# extract chromosomes 1-22 into separate files
for i in {1..22}
do
  echo "Creating vcf for chr${i}..."
  vcf="${out_prefix}_chr${i}.vcf.gz"
  bcftools view -r "chr${i}" -Oz "${in_vcf}" > "${vcf}"
  echo "${vcf} created."
  echo "Indexing VCF..."
  bcftools index "${vcf}"
  echo ""
done

# extract sex chromosomes (X and Y)
for i in "X" "Y"
do 
    echo "Creating vcf for chr${i}..."
    vcf="${out_prefix}_chr${i}.vcf.gz"
    bcftools view -r "chr${i}" -Oz "${in_vcf}" > "${vcf}"
    echo "${vcf} created."
    echo "Indexing VCF..."
    bcftools index "${vcf}"
    echo ""
done

# completion message
echo "Done."
date 
echo ""