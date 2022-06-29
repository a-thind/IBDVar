#!/bin/bash

# Anisha Thind, 5Jun2022

# Intended use:
# ./s06_annotate_vep.sh &> s06_annotate_vep.log

# stop at runtime errors
set -e

# start message
echo $0
date
echo ""

# files and folder
base_dir="/home/share"
data_dir="${base_dir}/data/s04_annotate"
basename="IHCAPX8_dragen_joint.clinvar.reheaded"
input_vcf="${data_dir}/${basename}.vcf.gz"
output_vcf="${data_dir}/${basename}.VEP.vcf.gz"

# Vep script
VEP="${base_dir}/tools/ensembl-vep/vep"
vep_report="${data_dir}/${basename}.VEP.html"
vep_cache_info="${data_dir}/${basename}.VEP.vep_cache_info.txt"
# VEP cache
vep_cache_dir="${base_dir}/tools/ensembl-vep/cache"
fasta="${vep_cache_dir}/homo_sapiens/106_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
cache_version="106"
cache_assembly="GRCh38"
# CADD 
cadd_dir="${base_dir}/resources/cadd"
cadd_snv="${cadd_dir}/whole_genome_SNVs.tsv.gz"
cadd_indels="${cadd_dir}/gnomad.genomes.r3.0.indel.tsv.gz"

echo "-- Files --"
echo ""
echo "Input VCF file: ${input_vcf}"
echo "Output VCF file: ${output_vcf}"
echo "VEP report: ${vep_report}"
echo ""
echo "-- VEP cache --"
echo ""
echo "VEP cache folder: ${vep_cache_dir}"
echo ""
echo "CADD annotation files:"
echo ""
echo "${cadd_snv}"
echo "${cadd_indels}"
echo ""

echo "Annotating with VEP..."
tabix -h -f "${input_vcf}" chr1 | "${VEP}" \
--format vcf \
--output_file "${output_vcf}" \
--stats_file "${vep_report}" \
--vcf \
--force_overwrite \
--compress_output bgzip \
--fork 4 \
--offline \
--cache \
--species homo_sapiens \
--dir "${vep_cache_dir}" \
--cache_version "${cache_version}" \
--assembly "${cache_assembly}" \
--pick \
--gencode_basic \
--sift b \
--polyphen b \
--symbol \
--max_af \
--plugin CADD,"${cadd_snv}","${cadd_indels}" 

#--fasta "${fasta}" 
#--hgvs, --numbers, --domains, --canonical, , --af, --af_1kg, --af_esp, --af_gnomad, --pubmed, --uniprot, --mane, --tsl, --appris, --gene_phenotype, 
#--everything, --check_existing, --total_length, --nearest symbol, --exclude_null_alleles, --check_ref, --regulatory, --variant_class
#--biotype \
#--mirna \ --ccds, --protein \

echo ""

# Index annotated VCF
echo "Indexing annotated VCF..."
bcftools index "${output_vcf}"
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