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

# If no threads specified or non-numeric then default is 4
if [[ -z "${threads}" || "${threads}" =~ !^[0-9]+ ]]; then
  threads=4
fi

# Filtering parameters 
if [[ -z "${GQ}" || "${GQ}" =~ !^[0-9]+ ]]; then
  echo "GQ parameter unspecified in config file... using default value (20)"
  GQ=20
fi

if [[ -z "${DP}" || "${DP}" =~ !^[0-9]+ ]]; then
  echo "DP parameter unspecified in config file... using default value (10)"
  DP=10
fi

if [ -z "${plink}" ]; then
   echo "Error: Missing PLINK path."
   exit 1
elif [ ! -d "${plink}" ]; then 
    echo "Error: Plink path is not a directory."
    exit 1
fi

if [ -z "${clinvar}" ]; then
   echo "Error: Missing ClinVar VCF file path."
   exit 1
elif [ ! -e "${clinvar}" ]; then
    echo "Error: ClinVar VCF path does not exist."
    exit 1
fi

# IBD detection parameters

if [ -z "${ibis}" ]; then
  echo "Error: Missing IBIS argument."
  exit 1
elif [ ! -d "${ibis}" ]; then
  echo "Error: IBIS path: ${ibis} is not a directory."
  exit 1  
fi

if [[ "${ibis_mt}" =~ !^[0-9]+ ]]; then
  echo "Error: 'ibis_mt' argument provided is non-numerical."
  exit 1
fi

# Check CADD folder path
if [ -z "${cadd}" ]; then
   echo "Error: Missing CADD path."
   exit 1
elif [ ! -d "${cadd}" ]; then 
    echo "Error: CADD path is not a directory."
    exit 1
fi

# check parameters
if [ ! -z "${genes}" ]; then
  if [ ! -e "${genes}" ]; then
      echo "Error: genes file specified ${genes} is not found."
      exit 1
  elif [[ "${genes}" =~ (*\.xlsx|*\.csv)  ]]; then
      echo "Error: genes of interest file is not an excel (.xlsx) or csv file."
      exit 1
  fi
fi