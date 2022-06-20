#!/bin/bash

# Anisha Thind, 16Jun2022

# stop script if non-zero exit status 
set -e
# stop script if variable value is unset
set -u
# stop entire pipe if non-zero status encountered
set -o pipefail

# starting message
echo $0
date
echo ""

# files and folders
base_dir="/home/share"
data_dir="${base_dir}/data/s07_select_haploblocks/truffle"
pheno_input="${data_dir}/pheno_input.txt"
pheno_colour="${data_dir}/pheno_colour.txt"
pheno_output="${data_dir}/ideogram"
ibd_file="${data_dir}/IHCAPX8_vcf_truffle.segments"
phenogram="${base_dir}/tools/phenogram/pheno_gram.rb"
human_genome="${phenogram%pheno_gram*}human_genome.txt"

# progress report
printf "PhenoGram %s\n" ` ruby $phenogram --version `
date
echo ""
echo "IBD segment file: ${ibd_file}"
echo "Phenogram input file: ${pheno_input}"
echo "Output ideogram file: ${pheno_output}.png"
echo ""


# create file for phenogram
echo "Creating ideogram showing IBD segments using PhenoGram..."
echo ""
echo "Creating input file for PhenoGram..."
echo ""

# check if IBD segment file is from TRUFFLE (.segments) or IBIS (.seg)
if [ "${ibd_file#*.}" == "segments" ]
then
    # extract required fields from IBD segments file
    awk 'BEGIN{ print "CHR", "\t", "POS", "\t", "END", "\t", "PHENOTYPE"} NR>1{ gsub("chr", "", $4); printf("%s\t%s\t%s\t%s-%s\n", $4,$5,$6,$2,$3) }' "${ibd_file}" > "${pheno_input}"
    echo "hello"
else
    awk 'BEGIN{ print "CHR", "\t", "POS", "\t", "END", "\t", "PHENOTYPE"} NR>1{ gsub("chr", "", $3); printf("%s\t%s\t%s\t%s-%s\n", $3,$4,$5,$1,$2) }' "${ibd_file}" > "${pheno_input}"
fi

echo "Phenogram input file created."
# assign colour number
awk '
BEGIN {
    printf("%s\t%s\t%s\t%s\t%s\n", "CHR", "POS", "END", "PHENOTYPE", "POSCOLOR")
}
NR>1{  
    if (!($4 in pheno)) {
        count++
        pheno[$4] = count
        } 
    $5 = pheno[$4]
    printf("%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5)
    }' "${pheno_input}" > "${pheno_colour}"
echo ""

# create ideogram
echo "Creating ideogram..."
ruby "${phenogram}" -i "${pheno_colour}" \
    -g "${human_genome}" \
    -c "exhaustive" \
    -p "proximity" \
    -T \
    -o "${pheno_output}"
echo ""

# completion message
echo "Done."
date
echo ""