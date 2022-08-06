# IBDVar
A tool for prioritising identity-by-descent (IBD) variants in Whole Genome Sequencing (WGS) data from families with rare heritable diseases. IBDVar consists of a variant prioritisation pipeline command-line program and an intereactive Shiny dashboard for starting the pipeline and visualising output.

## Variant Priorisation Pipelines
IBDVar can prioritise both short variants and structural variants (SV) from multi-sample VCF files generated from the Illumina DRAGEN Joint Genotype pipeline. Both pipelines can be initiated from the command-line or inside the Shiny dashboard "Start pipeline" tab.

## Short Variants
### Usage
```
./short_variants.sh -C pipeline.config [-m in_vcf.md5sum ]
```
#### Options:
  - ```-C```: config file (ending with .config) containing all parameters to execute the pipeline (required)
  - ```-m```: md5sum file to perform and and md5sum check on the input VCF file specified in the config file
