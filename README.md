# IBDVar
A tool for prioritising identity-by-descent (IBD) variants in Whole Genome Sequencing (WGS) data from families with rare heritable diseases. IBDVar consists of a variant prioritisation pipeline command-line program and an intereactive Shiny dashboard for starting the pipeline and visualising output.

## Variant Priorisation Pipelines
IBDVar can prioritise both short variants and structural variants (SV) from multi-sample VCF files generated from the Illumina DRAGEN Germline Pipeline. Both prioritisation pipelines can be initiated from the command-line or inside the Shiny dashboard "Start pipeline" tab.

## Short Variants
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


### Usage
```
./short_variants.sh -c pipeline.config [-m in_vcf.md5sum ]
```
#### Options:
  - ```-c```: config file (ending with .config) containing all parameters to execute the pipeline (required)
  - ```-m```: md5sum file to perform and and md5sum check on the input VCF file specified in the config file
  - ```-h```: help message with usage details and options

## Structural Variants
### Usage
```
./short_variants.sh -c pipeline.config
```
#### Options:
  - ```-c```: config file (ending with .config) containing all parameters to execute the pipeline (required)
  - ```-h```: help message with usage details and options
  
