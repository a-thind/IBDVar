if (!require(shinyjs)) {
  install.packages(shinyjs)
}

if (!require(reshape2)) {
  install.packages(reshape2)
}

library(shiny)
library(jsonlite)
library(shinyjs)
library(shinydashboard)
library(shinyFiles)
library(shinyFeedback)
library(biomaRt)
library(ideogram)
library(DT)
library(readxl)
library(dplyr)
library(purrr)
library(ggplot2)

options(shiny.maxRequestSize=1000*1024^2)



#' Filter dataframe by a given variable
#'
#' @param var variable to be filtered
#' @param in_val input value from widget in the UI
#'
#' @return
#' @export
#'
#' @examples
filter_variables <- function(var, in_val){
  if (is.factor(var)) {
    var %in% in_val
  } else if (is.character(var)) {
    # this clause is for variables with 2 or more levels
    map(in_val, function(x){grepl(x, var)}) %>% reduce(~.x|.y, .init=0)
  } else {
    # in case neither return TRUE
    TRUE
  }
}

# function to create a config file for the short variants pipeline
make_short_config <- function(in_vcf, out_dir, GQ, DP, 
                              MAF, ibis_mt1, ibis_mt2, mind, geno, threads) {
  # read in file paths
  # create a dataframe with parameters
  config_dir="scripts/config"
  config_path=file.path(config_dir, "pipeline_config.config")
  params_df <- data.frame(
    params=c("in_vcf", "out_dir", "GQ", "DP", "MAF", "ibis_mt1", "ibis_mt2", 
             "mind", "geno"),
    vals=c(in_vcf, out_dir, GQ, DP, MAF, ibis_mt1, ibis_mt2, mind, geno)
  )
  tools <- read.delim(file.path(config_dir, "tools_resources.cf"), 
                      comment.char = "#", sep="=", header = F)
  colnames(tools) <- c("params", "vals")
  # Concatenate parameter and tools together
  config <- rbind(params_df, tools)
  #print(config)
  # write a text file with "=" separator = config file
  write.table(config, config_path, col.names = F, row.names = F, sep="=", 
              quote = F)
}

# creates an SV config file from input form
sv_config <- function(in_vcf, out_dir, sv_ibis_seg, genes=NULL, sv_threads, email){
  config_dir="scripts/config"
  config_path=file.path(config_dir, "pipeline_sv.config")
  # if genes list is not provided
  if (is.null(genes)){
    params_df <- data.frame(
      params=c("sv_vcf", "out_dir", "ibd_seg", "threads"),
      vals=c(in_vcf, out_dir, sv_ibis_seg, sv_threads)
    )
  } else {
    params_df <- data.frame(
      params=c("sv_vcf", "out_dir", "ibd_seg", "threads", "genes"),
      vals=c(in_vcf, out_dir, sv_ibis_seg, sv_threads, genes)
    )
  }
 
  tools <- read.delim(file.path(config_dir, "tools_resources.cf"), 
                      comment.char = "#", sep="=", header = F)
  colnames(tools) <- c("params", "vals")
  # Concatenate parameter and tools together
  config <- rbind(params_df, tools)
  #print(config)
  # write a text file with "=" separator = config file
  write.table(config, config_path, col.names = F, row.names = F, sep="=", 
              quote = F)
}

#' function to parse consequences and get unique values (levels)
#'
#' @param var subsetted dataframe by a variable e.g. df$var
#'
#' @return
#' @export
#'
#' @examples
parse_levels <- function(var) {
  # paste the strings together
  concat_levels <- paste0(unique(var), collapse=", ")
  # split the strings
  parsed_levels <- unlist(strsplit(concat_levels, ', '))
  uniq_levels <- unique(parsed_levels)
  return(uniq_levels)
}

# create link for gene symbols to NCBI gene db
#' 
#'
#' @param column gene column of dataframe
#'
#' @return URL endpoint for access gene entry in NCBI
#' @export
#'
#' @examples
ncbi_gene <- function(column){
  ifelse(
    !is.na(column),
    paste0('<a href="https://www.ncbi.nlm.nih.gov/gene?term=(human[Organism]) AND ',
           column, '[Gene Name]">', column,'</a>')
  )
}


#' Convert multi-gene list strings to NCBI gene API request
#'
#' @param genes a string of genes delimited with a comma (row in SV table)
#'
#' @return
#' @export
#'
#' @examples
process_multigenes <- function(genes) {
  # split genes using "." delimiter
  split_genes <- unlist(strsplit(genes, ", "))
  # convert each gene to NCBI API route
  ncbi_genes <- sapply(split_genes, ncbi_gene)
  paste0(ncbi_genes, collapse=",")
  
}
