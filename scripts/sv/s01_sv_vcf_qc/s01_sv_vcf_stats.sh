#!/bin/bash
# s01_check_sv_vcf.sh - check
# Anisha Thind, 23Jun2022

# Intended use:
# ./s01_sv_vcf_stats.sh &> s01_sv_vcf_stats.log

# starting message
printf "Script: s01_sv_vcf_stats.sh\n\n"
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s00_source_data/IHCAPX8/s02_structural_variants"
in_vcf="${data_dir}/IHCAPX8_SV_dragen_joint.sv.vcf.gz"
csv="${data_dir}"
stats_dir="${base_dir}/data/sv/s01_sv_vcf_qc/bcfstats"
mkdir -p "${stats_dir}"
basename=$( basename ${in_vcf} .sv.vcf.gz ) 
stats_file="${stats_dir}/${basename}.vchk"


# Progress report
bcftools --version
date
echo ""

echo "Input VCF: ${in_vcf}"
echo "Output folder: ${stats_dir}"
echo ""

#echo "Indexing VCF..."
#bcftools index -f "${in_vcf}"
#echo ""

#echo "Generating bcfstats file..."
#bcftools stats -s - "${in_vcf}" > "${stats_file}"
#echo ""

#echo "Creating bcfstats plots"
#plot-vcfstats -s -p "${stats_dir}" "${stats_file}"
#echo ""

echo "Variant Counts"
echo "------------------"
# extract SVTYPE info from VCF
bcftools query -f "%SVTYPE\n"  "${in_vcf}" \
    | sort 
    | uniq -c 
    | awk '{ printf("%s %s %s:\t%s\n", $2, $3, $4, $1) }'

bcftools view -H "${in_vcf}" | awk ' END{ printf("Number of SV:\t%s\n\n", NR) } '




# Completion message
echo "Done."
date
echo ""
