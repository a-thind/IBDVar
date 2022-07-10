#!/bin/bash
# Anisha Thind, 7Jul2022

# filter variants with MAF above threshold
bcftools view -q "${MAF}":minor "${in_vcf}" \
  --threads "${threads}" \
  -Oz \
  -o "${filtered_vcf}"

zgrep -v "^#" "${filtered_vcf}" \
  | awk 'END{printf("Number of variants after filtering: %s\n\n", NR)}'