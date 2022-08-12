library(jsonlite)
library(dplyr)

# TODO:
# PR and SR table

# function to create a config file for the short variants pipeline
make_short_config <- function(in_vcf, out_dir, GQ, DP, 
                              MAF, ibis_mt1, ibis_mt2) {
  # read in file paths
  # create a dataframe with parameters
  config_dir="scripts/config"
  config_path=file.path(config_dir, "pipeline_config.config")
  params_df <- data.frame(
    params=c("in_vcf", "out_dir", "GQ", "DP", "MAF", "ibis_mt1", "ibis_mt2"),
    vals=c(in_vcf, out_dir, GQ, DP, MAF, ibis_mt1, ibis_mt2)
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


#' Create bins for CADD Phred values
#'
#' @param cadd CADD column of a dataframe
#'
#' @return
#' @export
#'
#' @examples
cadd_factors <- function(cadd) {
  cadd_factors <- as.factor(
    ifelse(is.na(cadd), 'NA',
           ifelse(cadd < 10, "< 10", 
                  ifelse(cadd < 20, "10 <= score < 20",
                         ifelse(cadd < 30, '20 <= score < 30', 
                                ifelse(cadd < 40, '30 <= score < 40',
                                       ifelse(cadd < 50, '40 <= score <50',
                                              ifelse(cadd < 60, '50 <= score< 60',
                                                     ifelse(cadd < 70, '60 <= score < 70',
                                                            ifelse(cadd < 80, '70 <= score < 80', 
                                                                   ifelse(cadd < 90, '80 <= score < 90',
                                                                          ifelse(cadd < 100, '90 <= score < 100'
                                                                          ))))))))))
    )
  )
  return(cadd_factors)
}

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
    var %in% in_val
  } else {
    # in case neither return null
    return(NULL)
  }
}

# Define server logic
server <- function(input, output, session) {
#-------------------------------------------------------------------------------
# Start pipeline tab
#-------------------------------------------------------------------------------
  options(shiny.maxRequestSize=1000*1024^2)
  # input vcf field
  sh_vcf <- reactive({
    req(input$sh_vcf)
    filename <- input$sh_vcf$name
    # check vcf has right extension (.vcf.gz)
    if (substr(filename, nchar(filename)-6, 
               nchar(filename))==".vcf.gz") {
      in_vcf=tools::file_path_as_absolute(input$sh_vcf$datapath)
      print(in_vcf)
    } else {
      validate("Input file is not a compressed VCF file (.vcf.gz).")
    }
  })
  
  # volumes set home directory
  volumes <- c(home="/home")
  shinyDirChoose(input, "sh_outdir", roots=volumes)
  
  
  sh_out_dir <- reactive({
    req(input$sh_outdir)
    out_dir <- parseDirPath(roots = volumes, input$sh_outdir)
  })
  
  output$sh_outdir_txt <- renderText({
    req(input$sh_outdir)
    paste0("Selected output folder: ", sh_out_dir())
  })
  
  # Genotype Quality
  GQ <- reactive({
    req(input$GQ)
    validate("Field required*")
    GQ <- input$GQ
    
  })
  # Depth
  DP <- reactive({
    req(input$DP)
    validate("Field required*")
    DP <- input$DP
  })
  # Minor Allele Frequency
  MAF <- reactive({
    req(input$MAF)
    validate("Field required*")
    MAF <- input$MAF
  })
  
  # IBD1 for IBIS
  ibis_mt1 <- reactive({
    validate("Field required*")
    ibis_mt1 <- input$ibis_mt1
  })
  # IBD2 for IBIS
  ibis_mt2 <- reactive({
    req(input$ibis_mt2)
    validate("Field required*")
    ibis_mt2 <- input$ibis_mt2
  })
  
  email <- reactive({
    req(input$email)
    validate("Please provide a valid email address for pipeline notifications.")
  })
  
  
  
  observeEvent(input$sh_start,{
    # make config file
    make_short_config(in_vcf=sh_vcf(), out_dir=sh_out_dir(), MAF=MAF(), GQ=GQ(), 
                      DP=DP(), 
                      ibis_mt1 = ibis_mt1(), ibis_mt2=ibis_mt2())
    showNotification("Short variants prioritisation pipeline started.",
                     type="message"
    )
    Sys.sleep(3)
    system("cd scripts; ./short_variants.sh -C config/pipeline_config.config", 
           wait=FALSE)
    
  })
  
  # start pipeline
  sv_vcf <- reactive({
    req(input$vcf)
    ext=tools::file_ext(input$vcf$name)
    switch(
      vcf.gz=input$vcf$datapath,
      vcf=input$vcf$datapath,
      validate("Input file is not a VCF file.")
    )
  })
  
#-------------------------------------------------------------------------------  
# Short Variants tab
#-------------------------------------------------------------------------------
  
  # read short variants csv
  col_types <- list(CHROM='f', ID='c', REF='c', ALT='c', FILTER='f', ALLELE='f',
                    CONSEQUENCE='c', IMPACT='f', SYMBOL='c', GENE='c',
                    FEATURE_TYPE='f', FEATURE='c', BIOTYPE='f', EXON='c',
                    INTRON='c', POS='d', QUAL='d', AC='d', AF='d', AN='d',
                    DP='d', FS='d', MQ='d', MQRANKSUM='d', QD='d',
                    READPOSRANKSUM='d', SOR='d', CDNA_POSITION='d',
                    CDS_POSITION='d', PROTEIN_POSITION='d', DB='f',
                    END='d', DISTANCE='d', FLAGS='f', MAX_AF_POPS='d',
                    SOMATIC='d', MOTIF_NAME='c', MOTIF_POS='d', HIGH_INF_POS='d',
                    MOTIF_SCORE_CHANGE='d', TRANSCRIPTION_FACTORS='c',
                    CADD_PHRED='d', HGVSC='c', HGVSP='c', PUBMED='c',
                    POLYPHEN_CALL='f',SIFT_CALL='f', CADD_RAW='d', SOR='d',
                    CODONS='c', STRAND='f', SIFT_SCORE='d', POLYPHEN_SCORE='d',
                    MAX_AF='d', EXISTING_VARIATION='c', HGNC_ID='c',
                    CLIN_SIG='f', SYMBOL_SOURCE='f', PHENO='f', AMINO_ACIDS='c',
                    NEAREST='c', HGVS_OFFSET='d', CLNSIG="f"
  )
  
  short_data <- reactive({
    req(input$short_tsv)
    ext=tools::file_ext(input$short_tsv$name)
    switch(ext,
           tsv=vroom::vroom(input$short_tsv$datapath, delim="\t",
                            col_names=TRUE, col_types = col_types),
           txt=vroom::vroom(input$short_tsv$datapath, delim="\t",
                            col_names=TRUE, col_types = col_types),
           validate("Invalid file: Please upload a tsv/text file")
    )
  })
  
  ibd_data <- reactive({
    req(input$ibd_seg)
    ext=tools::file_ext(input$ibd_seg$name)
    switch(ext,
           seg=input$ibd_seg$datapath,
           validate("Invalid file: Please upload an IBIS IBD segment file (.seg)")
    )
  })
  
  
  # render ideogram
  output$ideogram_plot <- renderIdeogram({
    ideogram({ibd_data()})
  })
  
  # reactive expression for ibd region filtering
  ibd_filter <- reactive({
    if (!is.null(input$chosenRegion$chr)) {
      short_data() %>%
        filter(CHROM==paste0('chr', input$chosenRegion$chr) &
                 (POS >= input$chosenRegion$start) &
                 (POS <= input$chosenRegion$stop) &
                 filters())
    } else {
      short_data() %>% filter(filters())
    }
  })
  
  # make filters
  # get vector of filter variables
  filters <- reactive({
    filter_variables(cadd_bins(), input$cadd_filter) &
      filter_variables(short_data()$IMPACT, input$impact_filter) &
      filter_variables(short_data()$SIFT_CALL, input$sift_filter) &
      filter_variables(short_data()$POLYPHEN_CALL, input$polyphen_filter) &
      filter_variables(short_data()$CONSEQUENCE, input$csq_filter) &
      filter_variables(short_data()$CLNSIG, input$clnsig_filter)
  })
  
  cadd_bins <- reactive({
    cadd_factors(short_data()$CADD_PHRED)
  })
  
  # filters panel triggered upon variant data loading
  output$filters  <- renderUI(
    tagList(
      checkboxGroupInput("csq_filter", "Consequence", 
                         selected=parse_levels(short_data()$CONSEQUENCE),
                         choices=parse_levels(short_data()$CONSEQUENCE)),
      checkboxGroupInput("clnsig_filter", "Clinical Signicance (ClinVar)", 
                         choices=levels(short_data()$CLNSIG),
                         selected=levels(short_data()$CLNSIG)),
      checkboxGroupInput("impact_filter", "VEP Impact",
                         selected=levels(short_data()$IMPACT),
                         choices=levels(short_data()$IMPACT)),
      checkboxGroupInput("sift_filter", "SIFT",
                         selected=levels(short_data()$SIFT_CALL),
                         choices=levels(short_data()$SIFT_CALL)),
      checkboxGroupInput("polyphen_filter", "PolyPhen", 
                         choices=levels(short_data()$POLYPHEN_CALL),
                         selected=levels(short_data()$POLYPHEN_CALL)),
      checkboxGroupInput("cadd_filter", "CADD Score", 
                         selected=levels(cadd_bins()), 
                         choices=levels(cadd_bins()))
      
    ))
  
  # render short variants table
  output$short_tab <- renderDT({ 
    DT::datatable(
      ibd_filter() %>%
        # create link for gene symbols to NCBI gene db
        mutate(
          SYMBOL=ifelse(
            !is.na(SYMBOL),
            paste0('<a href="https://www.ncbi.nlm.nih.gov/gene?term=(human[Organism]) AND ',
                   SYMBOL, '[Gene Name]">', SYMBOL,'</a>'), SYMBOL)) %>%
        mutate(
          RS=ifelse(
            !is.na(RS), paste0('<a href="https://www.ncbi.nlm.nih.gov/snp/?term=', 
                               RS, '">', RS, '</a>'), RS)
        ) %>%
        select(ID, RS, HGVSC, HGVSP,
               SYMBOL, CONSEQUENCE, MAX_AF, IMPACT,
               CADD_PHRED, POLYPHEN_CALL,
               SIFT_CALL, CLNSIG), 
      escape=FALSE,
      rownames=FALSE) 
  })
  
  # download variants results handler
  output$download <- downloadHandler(
    filename=function(){
      paste0(
        tools::file_path_sans_ext(input$short_tsv$name), ".tsv")
    },
    content=function(file){
      write.table({ibd_filter()}, file)
    })
  
  sv_data <- reactive({
    req(input$sv_tsv)
    ext=tools::file_ext(input$sv_tsv$name)
    switch(ext,
           tsv=vroom::vroom(input$sv_tsv$datapath, delim="\t", col_names = TRUE),
           txt=vroom::vroom(input$sv_tsv$datapath, delim="\t", col_names = TRUE),
           validate("Invalid file: Please upload a tsv/text file")
    )
  })
  
  output$sv_table <- renderDT({
    DT::datatable(
      sv_data() %>% 
        select(CHROM, START, END, ID, SV_TYPE, SV_LENGTH, CI_POS, CI_END, CIGAR, 
               GENE, OVERLAP))
  })
}

