#!/bin/bash
# s07_highlight_IBD.sh - Draw ideograms from IBD segment files from IBIS / TRUFFLE, using phenogram
# Anisha Thind, 16Jun2022

# Intended use:
# ./s07_highlight_IBD.sh out_dir phenogram genome &> s07_highlight_IBD.log
# Parameters:
#   $1: (out_dir) output folder
#   $2: (phenogram) phenogram folder
#   $3: (genome) Human genome text file with 3 columns (id, size, centromere), 
#       (for example seehttp://visualization.ritchielab.org/downloads/human_genome.txt)

# stop script if non-zero exit status 
set -euo pipefail

# starting message
printf "Script:\ts07_highlight_IBD.sh\n\n"
date
echo ""

for tool in ibis truffle 
do
printf "Drawing IBD segments in ideogram detected using ${tool}\n\n"
# files and folders
out_dir="${1}/s07_select_haploblocks/${tool}"
phenogram="${2}/pheno_gram.rb"
genome="${3}"
# output files
pheno_input="${out_dir}/pheno_input.txt"
pheno_colour="${out_dir}/pheno_colour.txt"
pheno_output="${out_dir}/ideogram"
ibd_file=$( find "${out_dir}" -name *.seg* )

if [ ! -e "${out_dir}" ]; then
    echo "Error: could not find ${out_dir} data folder."
    exit 1
fi

if [ -z "${phenogram}" ]; then
    echo "Error: Missing phenogram path."
    exit 1
elif [ ! -e "${phenogram}" ]; then
    echo "Error: Phenogram path does not exist."
    exit 1
fi

if [ -z "${genome}" ]; then
    echo "Error: Missing phenogram human genome text file."
    exit 1
elif [ ! -e "${genome}" ]; then
    echo "Error: Phenogram human genome text file not found."
    exit 1
fi


# progress report
printf "PhenoGram %s\n" $( ruby "${phenogram}" --version )  
date
echo ""

echo "IBD segment file: ${ibd_file}"
echo "Phenogram input file: ${pheno_input}"
echo "Genome file: ${genome}"
echo "Output ideogram file: ${pheno_output}.png"
echo ""

# check IBD segment file has been created
if [ -z "${ibd_file}" ]; then
  echo "IBD segment file not found in output folder."
  exit 1
fi


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
    -g "${genome}" \
    -c "exhaustive" \
    -p "proximity" \
    -T \
    -o "${pheno_output}"
echo ""

done

# completion message
echo "Done."
date
echo ""