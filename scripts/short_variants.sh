#!/bin/bash
# short_variants.sh - Script for starting short variants prioritisation pipeline
# Anisha Thind, 15June2022

# stop at runtime errors
set -e
# stop pipeline if non-zero status error
set -o pipefail

# starting message
echo ""
echo "Short Variants Pipeline"
echo ""
date
echo ""


# Check input
# generate usage message
usage()
{
   printf "Usage: ./short_variants.sh [-c short_variants.config] [-m md5sum_file.md5sum]  \n\n" 
   printf "Options:\n"
   printf "\t-c,\tconfiguration file\n"
   printf "\t-m,\tmd5sum file for input VCF file\n"
   printf "\t-h,\tdisplay this help message\n"
   exit 0;
}
# display usage message if no args given
[ $# -eq 0 ] && usage

echo "----------------------------------- Settings ----------------------------------"
echo ""
# check parameter arguments
while getopts "c:m:h" arg; do
   case "${arg}" in
      h)
         usage
         exit 0
      ;;
      c)
         if [ -e "${OPTARG}" ]; then
            filename=` basename ${OPTARG} `
            if [ "${filename##*.}" == "config" ]; then
                echo "Configuration file: ${OPTARG}"
                echo ""
                . "${OPTARG}"
                printf 'Input VCF file:\t%s\n\n' "${IN_VCF}"
                printf 'Output folder:\t%s\n\n' "${OUT_DIR}"
                printf 'Number of threads:\t%s\n\n' "${THREADS}"
            else
                echo "File ${OPTARG} is not a config file."
                exit 1
            fi
         else
            echo "File ${OPTARG} is not found."
            exit 1   
         fi
      ;;
      m)
         if [ -e "${OPTARG}" ]; then
            filename=` basename "${OPTARG}"`
            if [ "${filename##*.}" == "md5sum" ]; then
               echo "--------------------------------- MD5SUM Check --------------------------------"
               echo ""
               echo "MD5SUM file: ${OPTARG}"
               echo ""
               # create folder for qc
               DATA_DIR="${OUT_DIR}/s01_short_vcf_qc"
               mkdir -p "${DATA_DIR}"
               LOG="${DATA_DIR}/s01_check_vcf.log"
               # run md5sum check
               . s01_short_vcf_qc/s01_check_vcf.sh "${IN_VCF}" "${OPTARG}" "${DATA_DIR}" &> "${LOG}"
               cat "${LOG}"
               echo ""
            else
               echo "File ${OPTARG} is not a md5sum file." 
               exit 1
            fi
         else
            echo "File ${OPTARG} is not found."
            exit 1
         fi
      ;;
      \?)
         echo "Invalid option: ${OPTARG}"
         usage
         exit 1
      ;;
   esac
done

# Check VCF input
if [ ! -e "${IN_VCF}" ]
then
   echo "File ${IN_VCF} not found."
fi


echo "--------------------------- Computing BCFtools stats --------------------------"
echo ""
LOG="s01_short_vcf_qc/s02_bcfstats.log"
s01_short_vcf_qc/s02_bcfstats.sh "${IN_VCF}" "${OUT_DIR}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${OUT_DIR}/s01_short_vcf_qc"
echo ""
echo ""

echo "--------------------------- Variant Pre-processing ----------------------------"
echo ""
echo "Filtering VCF..."
DATA_DIR="${OUT_DIR}/s02_retain_pass_filter_vars"
LOG="s01_retain_pass_filter_vars.log"
s02_retain_pass_filter_vars/s01_retain_pass_filter_vars.sh "${IN_VCF}" "${OUT_DIR}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

# plot stats
LOG="s02_check_vcf_stats.log"
s02_retain_pass_filter_vars/s02_check_vcf_stats.sh "${OUT_DIR}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

# multiallelic site parsing
echo "Parsing multiallelic sites..."
echo ""
DATA_DIR="${OUT_DIR}/s03_split_MA_sites"
LOG="s01_split_MA_sites.log"
s03_split_MA_sites/s01_split_MA_sites.sh "${OUT_DIR}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

echo "------------------------------ Variant Annotation -----------------------------"
echo ""
DATA_DIR="${OUT_DIR}/s04_annotate_vars"
LOG="s01_add_variant_ids.log"
echo "Adding variant IDs..."
s04_annotate_vars/s01_add_variant_ids.sh "${OUT_DIR}" "${THREADS}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

echo "Checking ClinVar VCF..."
LOG="s02_check_clinvar_vcf.log"
s04_annotate_vars/s02_check_clinvar_vcf.sh "${OUT_DIR}" "${CLINVAR}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

echo "Renaming Chromosomes in ClinVar VCF..."
LOG="s03_rename_clinvar_chrs.log"
s04_annotate_vars/s03_rename_clinvar_chrs.sh "${OUT_DIR}" "${CLINVAR}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

printf "Annotating variants using ClinVar...\n\n"
LOG="s04_annotate_clinvar.log"
s04_annotate_vars/s04_annotate_clinvar.sh "${OUT_DIR}" "${CLINVAR}" "${THREADS}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

printf "Retaining only standard chromosomes in annotated VCF...\n\n"
LOG="s05_retain_std_chrs.log"
s04_annotate_vars/s05_retain_std_chrs.sh "${OUT_DIR}" "${THREADS}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

echo "-------------------------------- IBD Detection --------------------------------"
echo ""

echo "Creating Plink dataset..."
echo ""
DATA_DIR="${OUT_DIR}/s07_select_haploblocks"
LOG="s01_make_plink_dataset.log"
s07_select_haploblocks/s01_make_plink_dataset.sh "${OUT_DIR}" "${PLINK}" "${THREADS}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

echo "Adding recombination distances to Plink dataset..."
echo ""
LOG="s02_add_recomb_dists.log"
s07_select_haploblocks/s02_add_recomb_dists.sh "${OUT_DIR}" "${IBIS}" "${GENETIC_MAP}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

printf "Excluding SNPs lacking recombination distances...\n\n"
LOG="s03_retain_snps_with_recomb_dists.log"
s07_select_haploblocks/s03_retain_snps_with_recomb_dists.sh "${OUT_DIR}" "${PLINK}" "${THREADS}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

printf "Adding recombination distances to filtered Plink dataset...\n\n"
LOG="s04_add_recomb_dists.log"
s07_select_haploblocks/s04_add_recomb_dists.sh "${OUT_DIR}" "${IBIS}" "${GENETIC_MAP}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

printf "Detecting IBD segments with IBIS...\n\n"
LOG="s05_select_haploblocks_ibis.log"
s07_select_haploblocks/s05_select_haploblocks_ibis.sh "${OUT_DIR}" "${IBIS}" "${IBIS_MT}" "${THREADS}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

printf "Detecting IBD segments with TRUFFLE...\n\n"
LOG="s06_select_haploblocks_truffle.log"
s07_select_haploblocks/s06_select_haploblocks_truffle.sh "${OUT_DIR}" "${TRUFFLE}" "${IBS1M}" "${IBS2M}" "${THREADS}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

printf "Drawing ideograms highlighting IBD segments...\n\n"
LOG="s07_highlight_IBD.log"
s07_select_haploblocks/s07_highlight_IBD.sh "${OUT_DIR}" "${PHENOGRAM}" "${GENOME}" &> "${LOG}"
cat "${LOG}"
mv "${LOG}" "${DATA_DIR}"
echo ""

# completion messages
echo "Pipeline completed."
date
echo ""