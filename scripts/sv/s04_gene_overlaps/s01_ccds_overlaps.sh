#!/bin/bash
# s02_CCDS_overlaps.sh - find overlaps with coding sequences from CCDS
# Anisha Thind, 3Jul2022

# starting message
echo "Script: s02_overlaps_gene_list.sh"
date
echo ""

# set files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/sv/s03_annotate_sv"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.pass.read_support_updated.vcf.gz"
ccds_dir="${base_dir}/resources/ccds"
ccds="${ccds_dir}/CCDS.current.txt"
ccds_bed="${ccds%.txt}.bed"
filtered_ccds="${ccds_dir}/ccds_filtered.bed"
sort_ccds="${ccds%%.txt}_sort.txt"
threads=4
out_dir="${base_dir}/data/sv/s04_gene_overlaps"
mkdir -p "${out_dir}"
genes="/home/share/resources/gene_list.txt"
genes_bed="${out_dir}/genes.bed"
all_overlaps="${out_dir}/all_overlaps.bed"
gene_overlaps="${out_dir}/gene_overlaps.bed"

# progress report
bedtools --version
date
echo ""

# create a bed file from CCDS and retain only "Public CCDS"
echo -e "Creating bed file from CCDS...\n"
awk 'BEGIN{NR>1; FS="\t"; OFS="\t"} 
    $8~/^[0-9]+$/{
        if ($6~/Public/) {
            print $1, $8, $9, $2, $3, $4, $5, $6, $7, $10, $11
        }
    }' "${ccds}" > "${ccds_bed}"

# run ccds filter script
Rscript s02_filter_ccds.R -f "${ccds_bed}" -o "${ccds_dir}"

# extract matching gene list lines from filtered_bed
echo -e "Finding overlaps...\n"
grep -w -f "${genes}" "${filtered_ccds}" > "${genes_bed}"

# sort genes bed file
sort -k1,1 -k2,2n -k3,3n --parallel "${threads}" "${genes_bed}" > "${genes_bed}"

# decompress vcf
gunzip -c "${in_vcf}"
bedtools intersect -a "${in_vcf%.gz}" \
    -b "${filtered_ccds}" \
    -wa \
    -wb > "${all_overlaps}"

awk 'END{printf("Number of SV overlaps with sequences from CCDS: %s\n\n", NR)}' "${all_overlaps}"

# genes of interest
bedtools intersect -a "${in_vcf%.gz}" \
    -b "${genes_bed}" \
    -wa \
    -wb > "${gene_overlaps}"

awk 'END{printf("Number of SV overlaps with genes of interest: %s\n\n", NR)}' "${gene_overlaps}"


echo "Done."
date
echo ""