#!/bin/bash
# s04_annotate_sv.sh - Annotate SV

# stop at runtime errors
set -euo pipefail

echo "Script: s04_annotate_sv.sh"
date
echo ""

# set files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/sv/s03_filter_imprecise_vars"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.pass.precise.vcf.gz"
ref_seq_dir="${base_dir}/resources/refseq"
ref_seq="${ref_seq_dir}/hg38.ncbiRefSeq.gtf.gz"
dec_ref_seq="${ref_seq%%.gz}"
sort_ref_seq="${dec_ref_seq%%.gtf}_sort.gtf"
exons="${ref_seq_dir}/exons.bed"
threads=4
out_dir="${base_dir}/data/sv/s04_annotate_SV"
mkdir -p "${out_dir}"
overlaps="${out_dir}/sv_overlaps.vcf"
ibd_bed="${data_dir}/ibd.bed"
ibd_sv="${data_dir}/ibd_sv.bed"

# progress report
bedtools --version
date
echo ""

echo "Input VCF file: ${in_vcf}"
echo "GTF annotation file: ${ref_seq}"
echo "Output exons bed file: ${exons}"
echo ""

# decompress file
echo -e "Decompressing RefSeq GTF file...\n"
gzip -dkf "${ref_seq}"

# sort GTF file by chromosome then start position to accelerate overlap detection
echo -e "Sorting GTF file by chromosome then start position...\n"
sort -k 1,1 -k2,2n "${dec_ref_seq}" --parallel "${threads}" > "${sort_ref_seq}"

echo "Chromosome order in GTF..."
cat "${sort_ref_seq}" | awk '{ print $1 }' | uniq | awk 'NR<25{ print $1 }'
echo ""

echo -e "Extract exons from GTF...\n"
awk '$3~"exon"{print $0}' "${sort_ref_seq}" > "${exons}"

# interset to find overlaps between SV and gene features
echo -e "Finding overlaps between SV and genes...\n"
bedtools intersect -a "${in_vcf}" \
    -b "${exons}" \
    -wa \
    -wb > "${overlaps}" 

# Total number of overlaps:
awk 'END{ printf("Total number of overlaps: %d\n\n", NR)}' "${overlaps}"



# completion message
echo "Done."
date
echo ""