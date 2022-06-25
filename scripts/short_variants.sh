#!/bin/bash
# short_variants.sh - Script for starting short variants prioritisation pipeline
# Anisha Thind, 15June2022

# stop at runtime errors
set -e
# stop pipeline if non-zero status error
set -o pipefail

# starting message
echo -e "Short Variants Pipeline\n"
echo ""
date
echo ""

# Check input
# generate usage message
usage()
{
   echo -e "Usage: ./short_variants.sh -c short_variants.config [-m md5sum_file.md5sum]\n" 
   echo -e "Options:\n"
   echo -e "\t-C,\tconfiguration file\n"
   echo -e "\t-m,\tmd5sum file for input VCF file\n"
   echo -e "\t-h,\tdisplay this help message\n"
   exit 0;
}
# display usage message if no args given
[ $# -eq 0 ] && usage

echo -e "----------------------------------- Settings ----------------------------------\n"
# check parameter arguments
while getopts "C:m:h" arg; do
   case "${arg}" in
      h)
         usage
         exit 0
      ;;
      C)
         if [ -e "${OPTARG}" ]; then
            filename=$( basename ${OPTARG} ) 
            if [ "${filename##*.}" == "config" ]; then
                echo "Configuration file: ${OPTARG}"
                echo ""
                config="${OPTARG}"
                config/read_config.sh "${config}"
                . "${config}"
                printf 'Input VCF file:\t%s\n\n' "${in_vcf}"
                printf 'Output folder:\t%s\n\n' "${out_dir}"
                printf 'Number of threads:\t%s\n\n' "${threads}"
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
            filename=$( basename "${OPTARG}" ) 
            if [ "${filename##*.}" == "md5sum" ]; then
               md5sum="${OPTARG}"
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


echo -e "=================================== Pipeline ===================================\n"
# make log directory
log_dir="${out_dir}/logs"
mkdir -p "${log_dir}"
pipeline_log="${log_dir}/pipeline.log"

s01_short_vcf_qc/s00_start_qc.sh "${in_vcf}" "${out_dir}" \
   "${md5sum}" "${log_dir}" |& tee -a "${pipeline_log}"

#--------------------------- Variant Pre-processing ----------------------------
s02_retain_pass_filter_vars/s00_start_pre-processing.sh "${in_vcf}" "${out_dir}" \
   "${log_dir}" |& tee -a "${pipeline_log}"

# echo "------------------------------ Variant Annotation -----------------------------"
# echo ""
# DATA_DIR="${out_dir}/s04_annotate_vars"
# LOG="s01_add_variant_ids.log"
# echo "Adding variant IDs..."
# s04_annotate_vars/s01_add_variant_ids.sh "${out_dir}" "${threads}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# echo "Checking ClinVar VCF..."
# LOG="s02_check_clinvar_vcf.log"
# s04_annotate_vars/s02_check_clinvar_vcf.sh "${out_dir}" "${CLINVAR}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# echo "Renaming Chromosomes in ClinVar VCF..."
# LOG="s03_rename_clinvar_chrs.log"
# s04_annotate_vars/s03_rename_clinvar_chrs.sh "${out_dir}" "${CLINVAR}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# printf "Annotating variants using ClinVar...\n\n"
# LOG="s04_annotate_clinvar.log"
# s04_annotate_vars/s04_annotate_clinvar.sh "${out_dir}" "${CLINVAR}" "${threads}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# printf "Retaining only standard chromosomes in annotated VCF...\n\n"
# LOG="s05_retain_std_chrs.log"
# s04_annotate_vars/s05_retain_std_chrs.sh "${out_dir}" "${threads}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# echo "-------------------------------- IBD Detection --------------------------------"
# echo ""

# echo "Creating Plink dataset..."
# echo ""
# DATA_DIR="${out_dir}/s07_select_haploblocks"
# LOG="s01_make_plink_dataset.log"
# s07_select_haploblocks/s01_make_plink_dataset.sh "${out_dir}" "${PLINK}" "${threads}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# echo "Adding recombination distances to Plink dataset..."
# echo ""
# LOG="s02_add_recomb_dists.log"
# s07_select_haploblocks/s02_add_recomb_dists.sh "${out_dir}" "${IBIS}" "${GENETIC_MAP}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# printf "Excluding SNPs lacking recombination distances...\n\n"
# LOG="s03_retain_snps_with_recomb_dists.log"
# s07_select_haploblocks/s03_retain_snps_with_recomb_dists.sh "${out_dir}" "${PLINK}" "${threads}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# printf "Adding recombination distances to filtered Plink dataset...\n\n"
# LOG="s04_add_recomb_dists.log"
# s07_select_haploblocks/s04_add_recomb_dists.sh "${out_dir}" "${IBIS}" "${GENETIC_MAP}" &> "${LOG}"
# cat "${LOG}"
# mv "${LOG}" "${DATA_DIR}"
# echo ""

# printf "Detecting IBD segments with IBIS...\n\n"
# #LOG="s05_select_haploblocks_ibis.log"
# #s07_select_haploblocks/s05_select_haploblocks_ibis.sh "${out_dir}" "${IBIS}" "${IBIS_MT}" "${threads}" &> "${LOG}"
# #cat "${LOG}"
# #mv "${LOG}" "${DATA_DIR}"
# echo ""

# printf "Detecting IBD segments with TRUFFLE...\n\n"
# #LOG="s06_select_haploblocks_truffle.log"
# #s07_select_haploblocks/s06_select_haploblocks_truffle.sh "${out_dir}" "${TRUFFLE}" "${IBS1M}" "${IBS2M}" "${threads}" &> "${LOG}"
# #cat "${LOG}"
# #mv "${LOG}" "${DATA_DIR}"
# echo ""

# printf "Drawing ideograms highlighting IBD segments...\n\n"
# #LOG="s07_highlight_IBD.log"
# #s07_select_haploblocks/s07_highlight_IBD.sh "${out_dir}" "${PHENOGRAM}" "${GENOME}" &> "${LOG}"
# #cat "${LOG}"
# #mv "${LOG}" "${DATA_DIR}"
# echo ""

# completion messages
echo "Pipeline completed."
date
echo ""