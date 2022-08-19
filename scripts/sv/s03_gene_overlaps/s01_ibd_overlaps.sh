#!/bin/bash
# s01_ibd_overlaps.sh - find overlaps with coding sequences from CCDS and then genes of interest
# Anisha Thind, 3Jul2022

# starting message
echo "Script: s01_ibd_overlaps.sh"
date
echo ""

# set files and folders
out_dir="${1}"
data_dir="${out_dir}/s02_filter_sv"
in_vcf=$( find "${data_dir}" -name *.pass.vcf.gz )

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

# IBD files
ibd_seg="${4}"
ibd_bed="${out_dir}/ibd.bed"
ibd_sv_ccds="${out_dir}/ibd_filtered.txt"
ibd_var_ids="${out_dir}/ibd_variant_ids.txt"

# overlaps with CCDS (for annotation)
ibd_overlaps="${out_dir}/ibd_overlaps.bed"
# overlaps with genes
annot_ccds="${out_dir}/sv_ccds.bed"

if [ -z "${in_vcf}" ]; then
    "Error: missing pass filtered vcf file."
    exit 1
fi

# progress report
bedtools --version
date
echo ""

echo ""
echo "Input VCF: ${in_vcf}"
echo "Output Bed file: ${ibd_overlaps}"

# create bed file from segment file
awk -F"\t" '{OFS="\t"} NR>1{ print $3, $4, $5}' "${ibd_seg}" \
    | sed "s/^chr//" > "${ibd_bed}"

# decompress vcf
gunzip -c "${in_vcf}" > "${in_vcf%.gz}"
in_vcf="${in_vcf%.gz}"
echo -e "Checking chromosome notation in vcf...\n"
sed -i "s/^chr//" "${in_vcf}"


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

# report all genes with ccds annotation if there are overlaps
bedtools intersect \
    -a "${in_vcf}" \
    -b "${filtered_ccds}" \
    -filenames -wao > "${annot_ccds}"

awk 'BEGIN{FS="\t"; count=0} {
    if ($NF > 0) {count++}
    } 
    END{
        printf("\nNumber of SV overlaps with sequences from CCDS: %s\n\n", count)
    }' "${annot_ccds}"


# find which variants overlap with ibd regions
bedtools intersect \
    -a "${in_vcf}" \
    -b "${ibd_bed}" \
    -wa -wb > "${ibd_overlaps}"

zgrep -v "^#" "${ibd_overlaps}" \
    | awk 'END{printf("Number of SV overlaps with IBD regions: %s\n\n", NR)}'
# parse variant IDs
awk 'BEGIN{FS="\t"} {printf("%s\n", $3)}' "${ibd_overlaps}" > "${ibd_var_ids}"
echo ""
# Annotate IBD variants with CCDS
echo "${in_vcf}"
Rscript sv/s03_gene_overlaps/s03_annotate_ibd_sv.R "${in_vcf}" "${out_dir}" "${annot_ccds}" "${ibd_overlaps}"

echo "Done."
date
echo ""