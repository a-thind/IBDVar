#!/bin/bash
# s01_gene_overlaps.sh - find overlaps with coding sequences from CCDS and then genes of interest
# Anisha Thind, 3Jul2022

# starting message
echo "Script: s01_gene_overlaps.sh"
date
echo ""

# set files and folders
out_dir="${1}"
data_dir="${out_dir}/s02_filter_sv"
in_vcf=$( find "${data_dir}" -name *.pass.read_support.vcf.gz )

# ccds 
ccds_dir="${2}"
ccds="${ccds_dir}/CCDS.current.txt"
ccds_bed="${ccds%.txt}.bed"
filtered_ccds="${ccds_dir}/ccds_filtered.bed"
sort_ccds="${ccds%%.txt}_sort.txt"

threads="${3}"

# output folder
out_dir="${out_dir}/s03_gene_overlaps"
mkdir -p "${out_dir}"
# gene list files
genes="${4}"
genes_base="${genes%%.*}"
genes_bed="${genes_base}.bed"
genes_bed_sort="${genes_base}_sort.bed"
# overlaps with CCDS (for annotation)
all_overlaps="${out_dir}/all_overlaps.bed"
# overlaps with genes
gene_overlaps="${out_dir}/gene_overlaps.bed"


if [ -z "${in_vcf}" ]; then
    "Error: missing read supported filtered vcf file."
    exit 1
fi

# progress report
bedtools --version
date
echo ""

echo ""
echo "Input VCF: ${in_vcf}"
echo "Gene list file: ${genes}"
echo "Output Bed file: ${gene_overlaps}"

# create a bed file from CCDS and retain only "Public CCDS"
echo -e "Creating bed file from CCDS...\n"
awk 'BEGIN{NR>1; FS="\t"; OFS="\t"} 
    $8~/^[0-9]+$/{
        if ($6~/Public/) {
            print $1, $8, $9, $2, $3, $4, $5, $6, $7, $10, $11
        }
    }' "${ccds}" > "${ccds_bed}"

# run ccds filter script
Rscript sv/s03_gene_overlaps/s02_filter_ccds.R -f "${ccds_bed}" -o "${ccds_dir}"

# check if file is excel spreadsheet
if [ "${genes}" == "*.xlsx" ]; then
    in2csv "${genes}" > "${genes##.}.csv"
    genes="${genes##.}.csv"
fi

# extract matching gene list lines from filtered_bed
echo -e "Finding overlaps...\n"
grep -w -f "${genes}" "${filtered_ccds}" > "${genes_bed}"

# sort genes bed file
sort -k1,1 -k2,2n -k3,3n --parallel "${threads}" "${genes_bed}" > "${genes_bed_sort}"

# decompress vcf
gunzip -c "${in_vcf}" > "${in_vcf%.gz}"

in_vcf="${in_vcf%.gz}"
echo -e "Checking chromosome notation in vcf...\n"
sed -i "s/^chr//" "${in_vcf}"

bedtools intersect -a "${in_vcf}" \
    -b "${filtered_ccds}" \
    -wa \
    -wb > "${all_overlaps}"

awk 'END{printf("Number of SV overlaps with sequences from CCDS: %s\n\n", NR)}' "${all_overlaps}"

# genes of interest
bedtools intersect -a "${in_vcf}" \
    -b "${genes_bed_sort}" \
    -wa \
    -wb > "${gene_overlaps}"

awk 'END{printf("Number of SV overlaps with genes of interest: %s\n\n", NR)}' "${gene_overlaps}"


echo "Done."
date
echo ""