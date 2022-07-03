#!/bin/bash
# s03_explore_overlaps.sh - explore overlaps of SVs with genes

# start message
echo "Script: s03_explore_overlaps.sh"
date
echo ""

# set files and folders
base_dir="/home/share/"
data_dir="${base_dir}/data/sv/s04_annotate_sv"
overlaps="${data_dir}/sv_overlaps.vcf"
ibd="${base_dir}/data/output/IHCAPX8/s07_select_haploblocks/ibis/ibis.seg"
ibd_bed="${data_dir}/ibd.bed"
ibd_sv="${data_dir}/ibd_sv.bed"


echo "Input VCF file: ${overlaps}"
echo "Output IBD overlaps: ${ibd_sv}"
echo ""

awk -F"\t" '{print $17}' "${overlaps}" \
    | sort \
    | uniq -c \
    | awk 'END{ printf("Number of genes with SV overlaps: %s\n", NR) }'
echo ""

# Filter with IBD segments
echo -e "Creating bed file from IBD segments...\n"
awk 'NR>1{printf("chr%s\t%s\t%s\n",$3, $4, $5)}' "${ibd}" > "${ibd_bed}"

sed -i "s/^chr//g" "${ibd_bed}"

# check for overlaps with IBD regions
bedtools window -a "${overlaps}" -b "${ibd_bed}" > "${ibd_sv}"

awk 'END{ printf("Number of overlaps in IBD segments: %d\n", NR) }' "${ibd_sv}"
echo ""

echo -e "Number of SVs overlapping with genes in IBD regions:\n"
awk '{print $17}' "${ibd_sv}" \
    | sort \
    | uniq -c \
    | wc -l

echo ""

echo "Done."
date
echo ""
