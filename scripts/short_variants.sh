#!/bin/bash
# short_variants.sh - Main script for starting short variants prioritisation pipeline
# Anisha Thind, 15June2022

# stop at runtime errors and if pipe contains a non-zero status
set -eo pipefail

# function to generate usage message
usage()
{
   echo -e "Usage: ./short_variants.sh -c short_variants.config [-m md5sum_file.md5sum]\n" 
   echo -e "Options:\n"
   echo -e "  -c, configuration file\n"
   echo -e "  -m, md5sum file for input VCF file\n"
   echo -e "  -h, display this help message\n"
   exit 0;
}
# display usage message if no args given
[ $# -eq 0 ] && usage

# check parameter arguments
while getopts "c:m:h" arg; do
   case "${arg}" in
      h)
         usage
         exit 0
      ;;
      c)
         if [ -e "${OPTARG}" ]; then
            filename=$( basename ${OPTARG} ) 
            if [ "${filename##*.}" == "config" ]; then
               config="${OPTARG}"
               # reads config file
               utils/read_sh_config.sh "${config}"
               . "${config}"
               # make log directory
               log_dir="${out_dir}/logs"
               mkdir -p "${log_dir}"
               pipeline_log="${log_dir}/pipeline.log"
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

# shift processed params so next param is $1
shift $((OPTIND-1))

# if there is a remaining param display usage
if [ ! -z "${1}" ]; then
   usage
fi

# email pipeline started
if [ ! -z "${email}" ]; then
   mail -s "Short variants pipeline" "${email}" <<< "The Short variants prioritisation pipeline job has started."
fi

# starting message
echo -e "Short Variants Pipeline\n" &> "${pipeline_log}"
date &>> "${pipeline_log}"
echo -e "\n----------------------------------- Settings ----------------------------------\n" \
   &>> "${pipeline_log}"
echo -e "Configuration file: ${config}\n" \
   &>> "${pipeline_log}"
printf 'Input VCF file:\t%s\n\n' "${in_vcf}" \
   &>> "${pipeline_log}"
printf 'Output folder:\t%s\n\n' "${out_dir}" \
   &>> "${pipeline_log}"
printf 'Genotype quality threshold (GQ) per sample: %s\n\n' "${GQ}" \
   &>> "${pipeline_log}"
printf 'Read depth (FORMAT) threshold (GQ) per sample: %s\n\n' "${DP}" \
   &>> "${pipeline_log}"
printf 'PLINK "mind" parameter: %s\n\n' "${mind}" \
   &>> "${pipeline_log}"
printf 'PLINK "geno" parameter: %s\n\n' "${geno}" \
   &>> "${pipeline_log}"
printf 'PLINK "maf" parameter: %s\n\n' "${MAF}" \
   &>> "${pipeline_log}"
printf 'PLINK "maf" parameter: %s\n\n' "${MAF}" \
   &>> "${pipeline_log}"
printf 'IBIS "mt" parameter: %s\n\n' "${ibis_mt1}" \
   &>> "${pipeline_log}"
printf 'IBIS "mt2" parameter: %s\n\n' "${ibis_mt2}" \
   &>> "${pipeline_log}"
printf 'Genetic map: %s\n\n' "${genetic_map}" \
   &>> "${pipeline_log}"

if [ ! -z "${genes}" ]; then
   printf 'Genes list: %s\n\n' "${genes}" \
   &>> "${pipeline_log}"
fi

printf 'Number of threads: %s\n\n' "${threads}" \
   &>> "${pipeline_log}"
echo -e "Log files are created in ${log_dir}\n\n" \
   &>> "${pipeline_log}"

echo -e "=================================== Pipeline ===================================\n" &>> "${pipeline_log}"
short_variants/s01_short_vcf_qc/s00_start_qc.sh "${in_vcf}" \
    "${out_dir}" \
    "${md5sum}" \
    "${log_dir}" \
    &>> "${pipeline_log}"

#--------------------------- Variant Pre-processing ----------------------------
short_variants/s02_filter_short_vars/s00_start_pre-processing.sh "${in_vcf}" \
    "${out_dir}" \
    "${GQ}" \
    "${DP}" \
    "${threads}" \
    "${log_dir}" \
    &>> "${pipeline_log}"

#-------------------------------- IBD Detection --------------------------------
short_variants/s04_select_ibd_variants/s00_start_IBD_detection.sh "${out_dir}" \
   "${plink}" \
   "${mind}" \
   "${geno}" \
   "${MAF}" \
   "${threads}" \
   "${ibis}" \
   "${genetic_map}" \
   "${ibis_mt1}" \
   "${ibis_mt2}" \
   "${log_dir}" \
   &>> "${pipeline_log}" 

#------------------------------ Variant Annotation -----------------------------
short_variants/s05_annotate_vars/s00_start_annotation.sh "${out_dir}" \
    "${clinvar}" \
    "${vep}" \
    "${cadd}" \
    "${threads}" \
    "${log_dir}" \
    &>> "${pipeline_log}"

#------------------------------ Variant selection ------------------------------
# if a gene is provided then pass the argument to select variants
short_variants/s06_select_variants/s00_select_variants.sh "${out_dir}" \
"${max_af}" \
"${log_dir}" \
&>> "${pipeline_log}"

# if the user provides genes of interest list
if [ ! -z "${genes}" ]; then
   vars_out="${out_dir}/s06_select_variants/filtered_short_vars.tsv"
   genes_out="${vars_out%.tsv}_genes.tsv"
   # check if file is excel spreadsheet
   if [[ "${genes}" == *.xlsx ]]; then
      in2csv "${genes}" > "${genes%.*}.csv"
      genes="${genes%.*}.csv"
   fi

   echo -e "Filtering variants by genes of interest...\n" &>> "${pipeline_log}"
   awk 'NR==1{ print }' "${vars_out}" > "${genes_out}"
   if [ ! -z $( grep -w -f "${genes}" "${vars_out}" ) ]; then
      grep -w -f "${genes}" "${vars_out}" &>> "${genes_out}"
   fi
   
   awk 'END{printf("Number of variants overlapping with genes of interest: %s\n\n", NR-1)}' \
   "${genes_out}" &>> "${pipeline_log}"
fi

echo -e "--------------------------- Final Output ----------------------------\n" \
   &>> "${pipeline_log}"
final_dir="${out_dir}/final_output"
mkdir -p "${final_dir}"

ibd_seg="${out_dir}/s04_select_ibd_variants/ibis/ibis.seg"
# if gene list is provided put the gene filtered results as final output
if [ ! -z "${genes}" ]; then
   variants="${out_dir}/s06_select_variants/filtered_short_vars_genes.tsv"
else
   variants="${out_dir}/s06_select_variants/filtered_short_vars.tsv"
fi

cp "${variants}" "${final_dir}/"
cp "${ibd_seg}" "${final_dir}/"

echo "Final output is in ${final_dir}." &>> "${pipeline_log}"

# completion messages
echo "Pipeline completed." &>> "${pipeline_log}"
date &>> "${pipeline_log}"
echo "" &>> "${pipeline_log}"

if [ ! -z "${email}" ]; then
   echo -e "The SV prioritisation pipeline job is completed.\n
   To view the final output in the IBDVar application, 
   log into the server and type the following url in a web browser: 
   http://138.250.31.2:3737/anisha/IBDVar" |  mail -s "SV pipeline" "${email}"
fi