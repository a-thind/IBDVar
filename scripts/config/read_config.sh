#!/bin/bash
# read_config.sh - Read config file parameters

# Anisha Thind, 24Jun,2022

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
if [ -z "${in_vcf}" ]; then
   echo "Error: Missing input VCF file in config file."
   exit 1
elif [ ! -e "${in_vcf}" ]; then
    echo "Error: VCF file ${in_vcf} not found."
    exit 1
elif [[ "${in_vcf}" =~ (*\.vcf|*\.vcf\.gz)  ]]; then
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

if [ -z "${threads}" ]; then
    echo "Error: missing number of threads in config file."
    exit 1
fi
