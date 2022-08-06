#!/bin/bash
# s07_draw_ideogram.sh - draws ideogram using s07_draw_ideogram.R script
# Anisha Thind, 2July2022

# Example:
# ./s07_draw_ideogram.sh out_dir &> s07_draw_ideogram.sh
# Parameters:
#   $1: (outdir) output directory
#   $2: (scripts_dir) directory for IBD detection scripts - s07_select_haploblocks

# starting message
echo "Script: s07_draw_ideogram.sh"
date
echo ""

# files and folders
out_dir="${1}"
ibd_seg="${out_dir}/ibis/ibis.seg"
ideogram_dir="${out_dir}/ideogram"
ideogram="${ideogram_dir}/ideogram.html"
scripts_dir="${2}"

# check output dir
if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory argument."
   exit 1
elif [ ! -d "${out_dir}" ]; then
   echo "Error: output directory argument: ${out_dir} is not a directory."
   exit 1
fi
# check input IBD segment file
if [ ! -e "${ibd_seg}" ]; then
    echo "IBD (IBIS) segment file: '${ibd_seg}' not found."
    exit 1
fi

mkdir -p "${ideogram_dir}"

echo "Input IBD segment files: ${ibd_seg}"
echo "Output Ideogram: ${ideogram}"
echo "Output folder: ${ideogram_dir}"
echo ""

echo "Creating ideogram plot..."
Rscript "${scripts_dir}/s07_draw_ideogram.R" -f "${ibd_seg}" -o "${ideogram_dir}"
echo ""

echo -e "To view the interactive ideogram plot, open 'ideogram.html' in a web browser.\n"
# Completion message
echo "Done."
date
echo ""
