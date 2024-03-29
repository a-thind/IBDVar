# IBDVar
A tool for prioritising identity-by-descent (IBD) variants in Whole Genome Sequencing (WGS) data from families with rare heritable diseases. IBDVar consists of a variant prioritisation pipeline command-line program and an intereactive Shiny dashboard for starting the pipeline and visualising output.

## Table of Contents

- [Overview](https://github.com/a-thind/IBDVar#overview)
- [System Requirements](https://github.com/a-thind/IBDVar#system-requirements)
- [Variant Prioritisation Pipelines](https://github.com/a-thind/IBDVar#variant-priorisation-pipelines)
  - [Short Variants](https://github.com/a-thind/IBDVar#short-variants-command-line-start)
    - [Input VCF File](https://github.com/a-thind/IBDVar#input-vcf-file)
    - [Configuration Parameters](https://github.com/a-thind/IBDVar#configuration-parameters)
    - [Using GNU screen](https://github.com/a-thind/IBDVar#using-a-screen-to-run-the-short-variants-pipeline)
    - [Usage](https://github.com/a-thind/IBDVar#usage)
    - [Options](https://github.com/a-thind/IBDVar#options)
  - [Structural Variants](https://github.com/a-thind/IBDVar#structural-variants-command-line-start)
    - [Input VCF file](https://github.com/a-thind/IBDVar#input-vcf-file-1)
    - [Configuration Parameters](https://github.com/a-thind/IBDVar#configuration-file)
    - [Usage](https://github.com/a-thind/IBDVar#usage-1)
    - [Options](https://github.com/a-thind/IBDVar#options-1)
- [Shiny Dashboard](https://github.com/a-thind/IBDVar#shiny-dashboard)
  - [Start Pipeline](https://github.com/a-thind/IBDVar#start-pipeline)
  - [Short Variants](https://github.com/a-thind/IBDVar#short-variants)
  - [Structural Variants](https://github.com/a-thind/IBDVar#structural-variants)
  - [Configuration File](https://github.com/a-thind/IBDVar#configuration-file)
- [Questions, Feature Requests, Bug Reports and Issues](https://github.com/a-thind/IBDVar#questions-feature-requests-bug-reports-and-issues)
- [Licence](https://github.com/a-thind/IBDVar#licence)
- [Collaborators](https://github.com/a-thind/IBDVar#collaborators)

## Overview
The use of IBDVar follows a three step process:
<img src="https://user-images.githubusercontent.com/26285885/190384228-7eadde67-ebe1-474c-8812-5912973af9f3.jpg" width=800px/>
<img src="https://user-images.githubusercontent.com/26285885/190868992-cbbd1afa-6ca5-44ab-b307-e2315ee3fdad.png" width=800px>


The prioritisation pipeline is composed of two sub-pipelines (short variants and structural variants (SV)) that are started independently. Users can upload a multi-sample VCF file and configure the short variants or structural variants prioritisation pipeline in the Shiny dashboard or run the pipelines on a multi-sample VCF file at the command line using a configuration file. Once the pipeline has completed the output can be explored interactively in the corresponding pipeline tab in the Shiny dashboard. Unique to the tool, is the integration of IBD segment detection in variant prioritisation for WGS data. An overview of the key steps is shown below.<br/><br/>
<img src="https://user-images.githubusercontent.com/26285885/190386278-3a7186f8-1275-4e9f-9d0e-2bff54c8dd3a.jpg" width=800px/>


## System Requirements
For running the bash pipeline backend:
- Linux OS (developed and tested on Ubuntu 22 LTS)
- R (>=4.2)
- [BCFtools (1.15.1)](http://samtools.github.io/bcftools/)
- ClinVar VCF file (GRCh38)
- [IBIS (v1.20.9)](https://github.com/williamslab/ibis/)
- [Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)
- [CADD (v1.6) plugin resources (SNVs and indels)](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/)
- [CCDS](https://ftp.ncbi.nih.gov/pub/CCDS/current_human/CCDS.current.txt) (release number 22) text file

For the deploying the shiny dashboard, the following R dependencies are required:
- shiny
- shinydashboard
- shinyFiles
- shinyJS
- htmlwidgets
- dplyr
- jsonlite
- purrr
- readxl
- DT
- [ideogram](https://github.com/a-thind/ideogram)
- reshape2

To install these R packages, type the following in an R console:
```
install.packages(c("shiny", "shinydashboard", "shinyFiles", "shinyJS", "htmlwidgets", "dplyr", "jsonlite", "purrr", "readxl", "DT", "reshape2"))
```

To install the ideogram library, find the path of ideogram tarball file (.tar.gz) and type:
```
install.packages("path/to/ideogram_0.0.0.9000.tar.gz", type="source", repos=NULL)
```

## Variant Priorisation Pipelines
IBDVar can prioritise both short variants and structural variants (SV) from multi-sample VCF files generated from the Illumina DRAGEN Pipeline. Both prioritisation pipelines can be initiated from the command-line or inside the Shiny dashboard "Start pipeline" tab.

### Short Variants (Command line start)
#### Input VCF file
A multi-sample VCF file contained short variants (indels/ SNPs) called from the Illumina DRAGEN pipeline is used as input (see the Illumina website for details). The VCF file format should adhere to [version 4.2 specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). The pipeline expects chromosome naming to be prefixed with "chr" however, the tool checks for naming consistencies between the input VCF and the annotation resources implemented in the pipeline. 
#### Configuration Parameters
To run the short variants pipeline at the command line, you will need to create a configuration file with parameters (with "=" separating the parameter and its value) described in the table below:
<table>
<tbody>
<tr>
<td width="189">
<p><strong>Category</strong></p>
</td>
<td width="189">
<p><strong>Configuration parameter</strong></p>
</td>
<td width="483">
<p><strong>Description</strong></p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>General settings </strong></p>
</td>
<td width="189">
<p>in_vcf</p>
</td>
<td width="483">
<p>An input file path for the small variants VCF produced from Illumina DRAGEN Germline Pipeline.</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>out_dir</p>
</td>
<td width="483">
<p>An output directory path location to generate pipeline output</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>threads</p>
</td>
<td width="483">
<p>The number of threads (CPU) for executing the pipeline (default: 4)</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>QC filtering</strong></p>
</td>
<td width="189">
<p>GQ</p>
</td>
<td width="483">
<p>Minimum genotype quality threshold for each sample (default: 20)</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>DP</p>
</td>
<td width="483">
<p>Minimum (FORMAT) read depth threshold per sample (default:10)</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>MAF</p>
</td>
<td width="483">
<p>Minimum allele frequency for variants to be selected for the PLINK dataset</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>IBD detection</strong></p>
</td>
<td width="189">
<p>mind</p>
</td>
<td width="483">
<p>Maximum percentage of missing genotype data e.g., 0.1 excludes samples with &gt; 10% missing genotype data (default: 0.1)</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>geno</p>
</td>
<td width="483">
<p>Select variants with missing calling rates lower than the provided value (default: 0.1)</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>max_af</p>
</td>
<td width="483">
<p>Maximum allele frequency threshold for rare variants from the gnomAD, ESP or 100 genomes project populations. (Default: 0.05)</p>
<p>&nbsp;</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>ibis_mt1</p>
</td>
<td width="483">
<p>Minimum number of markers for IBIS to call a segment IBD1</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>ibis_mt2</p>
</td>
<td width="483">
<p>Minimum number of markers for IBIS to call a segment IBD2</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>genes</p>
</td>
<td width="483">
<p>A list of genes of interest for selecting variants in specified genes (optional)</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>Tools</strong></p>
</td>
<td width="189">
<p>tools_dir</p>
</td>
<td width="483">
<p>Optional tools base directory path for tools required by the pipeline</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>plink</p>
</td>
<td width="483">
<p>PLINK2 directory path</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>vep</p>
</td>
<td width="483">
<p>Vep executable file path</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>ibis</p>
</td>
<td width="483">
<p>Ibis directory path</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>Resources</strong></p>
</td>
<td width="189">
<p>resources</p>
</td>
<td width="483">
<p>Optional base directory path for resources</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>clinvar</p>
</td>
<td width="483">
<p>ClinVar VCF file path</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>genetic_map</p>
</td>
<td width="483">
<p>The file path for the genetic recombination map for the human genome</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>cadd</p>
</td>
<td width="483">
<p>CADD plugin resource directory path</p>
</td>
</tr>
</tbody>
</table>

Click [here](./scripts/config/example_short.config) for an example of a short variants config file.

### Using a screen to run the short variants pipeline
As the short variants pipeline can take a few hours to complete, it is highly recommended to run the pipeline in a Linux GNU screen to prevent abrupt termination of the pipeline, for example, in the event of a connection drop or a sudden SSH session termination. To install Linux GNU Screen on Ubuntu / Debian systems:
```
sudo apt update
sudo apt install screen 
```
On CentOS/Fedora type:
```
sudo yum install screen
```
To create a screen type ```screen``` in the terminal, or create a named screen by typing the following:
```
screen -S <screen_name>
```
Attach the screen to the terminal as follows:
```
screen -r <screen_name>
```
Once the screen is attached execute the pipeline as described in the [usage](#usage) section below.

After attaching the screen to the terminal and initiating the short variants pipeline, detach the screen by pressing ```CtrA``` and ```d```, or typing in the terminal:
```
screen -d <screen_name>
```
This will allow exiting of the terminal window without terminating the pipeline. To reattach the screen, simply type in the terminal:
```
screen -r <screen_name>
```


#### Usage
```
./short_variants.sh -c pipeline.config [-m in_vcf.md5sum ]
```
##### Options:
  - ```-c```: config file (ending with .config) containing all parameters to execute the pipeline (required)
  - ```-m```: md5sum file to perform and and md5sum check on the input VCF file specified in the config file
  - ```-h```: help message with usage details and options

### Structural Variants (Command line start)
#### Input VCF File
A multi-sample VCF file contained structural variants called (using Manta) from the Illumina DRAGEN pipeline is used as input (see the Illumina website for details). The VCF file format should adhere to [version 4.2 specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). The pipeline expects chromosome naming to be prefixed with "chr" however, the tool checks for naming consistencies between the input VCF and the annotation resources implemented in the pipeline. 
#### Configuration File
To start the structural variants pipeline at the command line, you will need to create a configuration file using the parameters specified in the table below:
<table width="860">
<tbody>
<tr>
<td width="189">
<p><strong>Category</strong></p>
</td>
<td width="189">
<p><strong>Configuration parameter</strong></p>
</td>
<td width="482">
<p><strong>Description</strong></p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>General settings</strong></p>
</td>
<td width="189">
<p>sv_vcf</p>
</td>
<td width="482">
<p>Input VCF file path</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>out_dir</p>
</td>
<td width="482">
<p>Directory path for pipeline output</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>threads</p>
</td>
<td width="482">
<p>Number of threads (CPU)</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>Variant selection</strong></p>
</td>
<td width="189">
<p>ibd_seg</p>
</td>
<td width="482">
<p>IBD segment file path (from the short variants pipeline) for selecting SV in IBD segments.</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>genes</p>
</td>
<td width="482">
<p>A list of genes of interest to be used to filter variants.</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>Tools</strong></p>
</td>
<td width="189">
<p>tools_dir</p>
</td>
<td width="482">
<p>(Optional) base directory for tools</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>Resources</strong></p>
</td>
<td width="189">
<p>resources</p>
</td>
<td width="482">
<p>The base directory for resources (optional)</p>
</td>
</tr>
<tr>
<td width="189">
<p><strong>&nbsp;</strong></p>
</td>
<td width="189">
<p>ccds</p>
</td>
<td width="482">
<p>CCDS directory path</p>
</td>
</tr>
</tbody>
</table>

Click [here](./scripts/config/example_sv.config) for an example of a structural variants config file.

#### Usage
```
./structural_variants.sh -c pipeline.config
```
##### Options:
  - ```-c```: config file (ending with .config) containing all parameters to execute the pipeline (required)
  - ```-h```: help message with usage details and options
  
# Shiny Dashboard
The shiny dashboard allows users to start prioritisation pipelines for short or strucutral variants and to analyse the output interactively.

To start the Shiny Dashboard in the Cranfield Univeristy server, log into the Linux server deploying the tool and type the application URL (can be requested from the author) in a web-browser. (Note that development and testing was performed using the Google Chrome browser so performance may vary with other browsers.) The shiny dashboard can also be started in RStudio however it is not recommended, since most of views have been configured for browser display and may affect performance of the tool.

The "Start Pipeline" tab will be open first by default. 
## Start Pipeline
In the "Start Pipeline" tab you can start the short variants or structural variants pipeline by selecting an input VCF file, output folder for results and configuring  parameters listed in the respective pipeline box. <br/><br/>
Once parameters have been specified, click ```Start``` in the respective pipeline box to run the pipeline. A notification message should appear in the bottom right corner indicating pipeline initiation.<br/><br/>
![start_pipeline](https://user-images.githubusercontent.com/26285885/188452229-47ff8f44-dff3-4866-8a0b-0fc82c39b740.png)

## Short Variants
In the "Short Variants" tab you can explore the short variants pipeline output interactively. <br/>
The tab features:
- a "Files" box to upload the following files which are located in the "final_output" folder of the output folder specified at run-time of the pipeline: 
  1. A <b>prioritised and annotated list of variants</b> produced from the short variants prioritisation pipeline.
  2. An IBIS <b>IBD segment file</b> produced from the pipeline
  - An optional file containing  <i>list of genes of interest</i> can also be uploaded to filter the variants by these genes

- <b>Interactive variants table</b> - users can filter, sort, search and download a TSV file of variants reported in the table.
- <b>Filters panel</b> - contains a series of checkboxes to filter variants by CADD score, predicted consequence, SIFT and PolyPhen calls, clinical significance (ClinVar) and VEP predicted impact (loss of function etc.)
- <b>Interactive ideogram</b> - filters variants in the interactive data table below by the IBD region clicked by the user. A tool-tip reporting the chromosome number, start and end position of a given IBD region is displayed when a user hovers over an IBD region.<br/><br/>
![short_variants_tab](https://user-images.githubusercontent.com/26285885/188453809-bf266487-f180-4519-8aad-8986626d5f25.png)
- "Summary" box summarising: 
  - total number of variants
  - number of pathogenic variants identified by ClinVar
  - number of detected IBD segments, the total number of deleterious missense variants predicted by SIFT, PolyPhen and CADD
  - number of loss of function variants

## Structural Variants
In the "Structural Variants" tab, the prioritised SV calls from the pipeline can be interactively explored using filters and an interactive data table.<br/><br/>
SV tab features include:
- "Files" box for uploading the prioritised list of SV calls (.tsv) file
- <b>Interactive table</b> of variants that can filtered, sorted, searched and downloaded as a TSV file.
- "Summary" tab providing summary statistics on the various counts of SV types and also the mean SV lengths.
- <b>Filters panel</b> containing checkboxes to filter the variants table by: SV type, chromosome number, precision of breakpoints of called SVs and genes of interest.
<br/><br/>
![sv_tab (1)](https://user-images.githubusercontent.com/26285885/188453867-61f3e77d-2b45-48bb-980a-d589daad8d41.png)

# Questions, Feature Requests, Bug Reports and Issues
For any questions, feature requests, bug reports or issues regarding the latest version of IBDVar, please click on the "[issues](https://github.com/a-thind/IBDVar/issues)" tab present at the top-left of the GitHub repository page.

# Licence
[MIT](LICENSE)

# Collaborators
This codebase was developed as part of an MSc thesis project (MSc Applied Bioinformatics, Cranfield University 2021-2022) under the supervision of Dr Alexey Larionov.
