library(jsonlite)
library(dplyr)

# function to parse consequences and get unique values (levels)
parse_levels <- function(df, var) {
  # paste the strings together
  concat_levels <- paste0(unique(df$var), collapse=", ")
  # split the strings
  parsed_levels <- unlist(strsplit(concat_levels, ', '))
  uniq_levels <- unique(parsed_levels)
  return(uniq_levels)
}

# function to create select inputs for factor variables
ui_filters <- function(var) {
  levels <- levels(var)
  selectInput(var, var, choices=levels, selected = levels, multiple = TRUE)
}

#' Filter dataframe by a given variable
#'
#' @param var variable to be filtered
#' @param val input value from widget in the UI
#'
#' @return
#' @export
#'
#' @examples
filter_variables <- function(var, val){
  if (is.factor(var)) {
    var %in% val
  } else if (is.numeric(var)) {
    # filter variable between slider limits
    var >= val[1] & var <= val[2] & !is.na(var)
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
  
  output$cadd_test <- renderText({
    input$cadd_score
  })
  
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
                    NEAREST='c', HGVS_OFFSET='d'
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
    #filter_variables(short_data()$CADD_PHRED, input$cadd_filter) &
    filter_variables(short_data()$IMPACT, input$impact_filter) &
      filter_variables(short_data()$SIFT_CALL, input$sift_filter) &
      filter_variables(short_data()$POLYPHEN_CALL, input$polyphen_filter)
  })
  
  output$filters  <- renderUI(
    tagList(
      tags$h3("Filters"),
      sliderInput("cadd_filter", "CADD Score", value=20, min=1, max=100),
      checkboxGroupInput("impact_filter", "VEP Impact",
                         selected=levels(short_data()$IMPACT),
                         choices=levels(short_data()$IMPACT)),
      checkboxGroupInput("sift_filter", "SIFT",
                         selected=levels(short_data()$SIFT_CALL),
                         choices=levels(short_data()$SIFT_CALL)),
      checkboxGroupInput("polyphen_filter", "PolyPhen", 
                         choices=levels(short_data()$POLYPHEN_CALL),
                         selected=levels(short_data()$POLYPHEN_CALL))
    )
  )
  
  output$txt <- renderText({
    impact <- paste(input$impact_filter, collapse = ", ")
    paste("You chose", impact)
  })
  
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
               SIFT_CALL, CLNSIG, everything()), 
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

