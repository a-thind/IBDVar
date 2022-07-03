#!/bin/bash
# s04_annotate_sv.sh - Annotate SV

# stop at runtime errors
set -euo pipefail

echo "Script: s04_annotate_sv.sh"
date
echo ""

# set files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/sv/s04_annotate_sv"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.pass.precise.vcf"
ccds_dir="${base_dir}/resources/ccds"
ccds="${ccds_dir}/CCDS.current.txt"
ccds_bed="${ccds%.txt}.bed"
filtered_ccds="${ccds_bed%.bed}_filtered.bed"
sort_ccds="${ccds%%.txt}_sort.txt"
threads=4
out_dir="${base_dir}/data/sv/s04_annotate_sv"
mkdir -p "${out_dir}"
overlaps="${out_dir}/sv_overlaps.vcf"


# progress report
bedtools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo "CCDS annotation file: ${ccds}"
echo "Output overlaps file: ${overlaps}"
echo ""


# sort GTF file by chromosome then start position to accelerate overlap detection
#echo -e "Sorting GTF file by chromosome then start position...\n"
#sort -k 1,1n -k8,8n "${dec_ref_seq}" --parallel "${threads}" > "${sort_ref_seq}"

# create a bed file from CCDS and retain only "Public CCDS"
echo "Creating bed file from CCDS..."
awk 'BEGIN{NR>1; FS="\t"; OFS="\t"} 
    $8~/^[0-9]+$/{
        if ($6~/Public/) {
            print $1, $8, $9, $2, $3, $4, $5, $6, $7, $10, $11
        }
    }' "${ccds}" > "${ccds_bed}"

# interset to find overlaps between SV and gene features
echo -e "Finding overlaps between SV and genes...\n"
bedtools intersect -a "${in_vcf}" \
    -b "${ccds_bed}" \
    -wa \
    -wb > "${overlaps}" 

# Total number of overlaps:
awk 'END{ printf("Total number of overlaps: %d\n\n", NR)}' "${overlaps}"

# completion message
echo "Done."
date
echo ""