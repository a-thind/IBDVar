library(jsonlite)
library(dplyr)

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
           ifelse(cadd < 20, "< 20",
              ifelse(cadd < 30, '< 30', 
                     ifelse(cadd < 40, '< 40',
                                ifelse(cadd < 50, '<50',
                                       ifelse(cadd < 60, '< 60',
                                              ifelse(cadd < 70, '< 70',
                                                     ifelse(cadd < 80, '< 80', 
                                                            ifelse(cadd < 90, '< 90',
                                                                   ifelse(cadd < 100, '< 100'
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
  options(shiny.maxRequestSize=1000*1024^2)
  
  shinyDirChoose(input, "sh_outdir", roots=c(wd='.'))
  
  # start pipeline
  svvcf <- reactive({
    req(input$vcf)
    ext=tools::file_ext(input$vcf$name)
    switch(
      vcf.gz=input$vcf$datapath,
      vcf=input$vcf$datapath,
      validate("Input file is not a VCF file.")
    )
  })
  
  # Short Variants tab
  #-----------------------
  
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
      checkboxGroupInput("cadd_filter", "CADD Score", 
                         selected=levels(cadd_bins()), 
                         choices=levels(cadd_bins())),
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
                         selected=levels(short_data()$POLYPHEN_CALL))
      
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
  
}

