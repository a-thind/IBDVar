#!/bin/bash
# s06_annotate_vep.sh - annotate short variants using VEP
# Anisha Thind, 5Jun2022

# Example:
# ./s06_annotate_vep.sh &> s06_annotate_vep.log
# Parameters:
#   $1: (out_dir) output directory
#   $2: (vep) vep executable path
#   $3: (cadd) CADD directory
#   $4: (threads) number of threads


# stop at runtime errors
set -euo pipefail

# start message
echo -e "Script: s06_annotate_vep.sh\n"
date
echo ""

# files and folder
out_dir="${1}"
in_vcf="/home/share/data/output/IHCAPX8/s06_select_haploblocks/ibis/filtered_ibd.sorted.vcf.gz"
vep="${2}"
cadd_dir="${3}"
threads="${4}"
echo "Outdir ${out_dir}"

# check output
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi

# check update VCF exists
if [ -z "${in_vcf}" ]; then
  echo "Updated ClinVar VCF file not found."
  exit 1
fi

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ !^[0-9]+ ]]; then
  threads=4
fi

basename=$( basename "${in_vcf}" .vcf.gz )
output_vcf="${out_dir}/${basename}.vep.vcf.gz"
vep_report="${out_dir}/${basename}.vep.html"
vep_cache_info="${out_dir}/${basename}.vep.vep_cache_info.txt"
cache_dir="/home/share/tools/ensembl-vep/cache"
# VEP cache
cache_version="106"
cache_assembly="GRCh38"
# CADD 
cadd_snv="${cadd_dir}/whole_genome_SNVs.tsv.gz"
cadd_indels="${cadd_dir}/gnomad.genomes.r3.0.indel.tsv.gz"


echo ""
echo "Input VCF file: ${in_vcf}"
echo "Output VCF file: ${output_vcf}"
echo "VEP report: ${vep_report}"
echo ""
echo "CADD annotation files:"
echo ""
echo "${cadd_snv}"
echo "${cadd_indels}"
echo ""

echo "Annotating with VEP..."
"${vep}" \
  -i "${in_vcf}" \
  --format vcf \
  --output_file "${output_vcf}" \
  --stats_file "${vep_report}" \
  --vcf \
  --force_overwrite \
  --compress_output bgzip \
  --fork "${threads}" \
  --offline \
  --cache \
  --dir "${cache_dir}" \
  --species homo_sapiens \
  --cache_version "${cache_version}" \
  --assembly "${cache_assembly}" \
  --buffer_size 50 \
  --pick \
  --gencode_basic \
  --sift b \
  --polyphen b \
  --symbol \
  --max_af \
  --hgvs \
  --pubmed \
  --check_existing \
  --total_length \
  --nearest symbol \
  --regulatory \
  --check_ref \
  --exclude_null_alleles \
  --plugin CADD,"${cadd_snv}","${cadd_indels}" 

#--fasta "${fasta}" 
#--hgvs, --numbers, --domains, --canonical, , --af, --af_1kg, --af_esp, --af_gnomad, --pubmed, --uniprot, --mane, --tsl, --appris, --gene_phenotype, 
#--everything, --check_existing, --total_length, --nearest symbol, --exclude_null_alleles, --check_ref, --regulatory, --variant_class
#--biotype \
#--mirna \ --ccds, --protein \

echo ""

# Index annotated VCF
echo "Indexing annotated VCF..."
bcftools index -f "${output_vcf}"
echo ""

# VEP annotations
echo "Included VEP annotations:"
echo ""
bcftools +split-vep "${output_vcf}" -l
echo ""

# completion message
echo "Done."
date
echo ""