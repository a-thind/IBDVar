#!/bin/bash
# read_sv_ config.sh - Read config file parameters

# Anisha Thind, 6Aug,2022

config="${1}"

if [ ! -e "${config}" ]; then
   echo "Config file not found."
   exit 1
   if [ "${config%%.config}" == .config ]; then
        echo "File ${config} is not a config file."
        exit 1
    fi
fi

source "${config}"

# check parameters
if [ -z "${sv_vcf}" ]; then
   echo "Error: Missing input VCF file in config file."
   exit 1
elif [ ! -e "${sv_vcf}" ]; then
    echo "Error: VCF file ${sv_vcf} not found."
    exit 1
elif [[ "${sv_vcf}" =~ (*\.vcf|*\.vcf\.gz)  ]]; then
    echo "Error: input file is not a VCF file."
    exit 1
fi

if [ -z "${out_dir}" ]; then
   echo "Error: Missing output directory path in config file."
   exit 1
elif [ ! -d "${out_dir}" ]; then
    echo "${out_dir} is not a directory."
    exit 1
fi

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ !^[0-9]+ ]]; then
  threads=4
fi

# Filtering parameters 
if [[ -z "${PR}" || "${PR}" =~ !^[0-9]+ ]]; then
  echo "PR parameter unspecified in config file... using default value (8)"
  PR=8
fi

if [[ -z "${SR}" || "${SR}" =~ !^[0-9]+ ]]; then
  echo "SR parameter unspecified in config file... using default value (0.15)"
  SR=0.15
fi

# check ccds directory path
if [ -z "${ccds}" ]; then
   echo "Error: Missing CCDS folder path."
   exit 1
elif [ ! -d "${ccds}" ]; then 
    echo "Error: CCDS path is not a directory."
    exit 1
fi

if [ -z "${ibd_seg}" ]; then
   echo "Error: Missing IBD segment file path."
   exit 1
elif [ ! -e "${ibd_seg}" ]; then
    echo "Error: IBD segment file path does not exist."
    exit 1
fi


if [ -z "${genes}" ]; then
   echo "Error: Missing list of genes of interest file path."
   exit 1
elif [ ! -e "${genes}" ]; then
    echo "Error: genes list file path does not exist."
    exit 1
fi
