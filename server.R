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



# Define server logic
server <- function(input, output, session) {
#-------------------------------------------------------------------------------
# Start pipeline tab
#-------------------------------------------------------------------------------
  # input vcf field
  volumes <- getVolumes()
  shinyFileChoose(input, "sh_vcf", roots=volumes())
  sh_vcf_name <- reactive({
    if (!is.null(input$sh_vcf)){
      infile <- parseFilePaths(roots=volumes(), input$sh_vcf)$datapath
      if (length(infile) > 0){
        # check vcf has right extension (.vcf.gz)
        ifelse(substr(infile, nchar(infile)-6,
                      nchar(infile))==".vcf.gz",
               infile,
               validate("Invalid file: Input file is not a compressed VCF file (.vcf.gz).")
        )
      }
    }
    
  })
  
  
  output$sh_vcf_name <- renderText({
    req(input$sh_vcf)
    if (!is.null(input$sh_vcf)) {
      infile <- paste("Input VCF file: ", 
                      parseFilePaths(roots=volumes(), input$sh_vcf)$name)
    } else {
      "Input VCF file: "
    }
    
  })
  
  output$sh_vcf_label <- renderText({
    "<b>Upload a compressed short variants VCF file (.vcf.gz)</b>"
  })

  shinyDirChoose(input, "sh_outdir", roots=volumes())
  
  
  sh_out_dir <- reactive({
    req(input$sh_outdir)
    out_dir <- substring(parseDirPath(roots=volumes(), input$sh_outdir), 2)
  })
  
  output$sh_outdir_txt <- renderText({
    req(input$sh_outdir)
    prefix <- "Selected output folder: "
    if (!is.null(input$sh_outdir)){
      paste0(prefix, sh_out_dir()) 
    } else {
      prefix
    }
  })
  
  output$sh_outdir_label <- renderText({
    "<b>Select an output folder</b>"
  })
  
  shinyFileChoose(input, "sh_start_genes", roots=volumes())
  
  sh_start_genes <- reactive({
    if (!is.null(input$sh_start_genes)){
      infile <- parseFilePaths(roots=volumes(), input$sh_start_genes)$datapath
      if (length(infile) > 0){
        ext=tools::file_ext(infile)
        # check segment file 
        switch(ext,
               xlsx=infile,
               txt=infile,
               validate("Invalid file: Input file is not an excel (.xlsx) or text file (.txt")
        )
      }
    }
  })
  
  output$sh_start_genes_name <- renderText({
    req(input$sh_start_genes)
    if (!is.null(input$sh_start_genes)) {
      infile <- paste("Filename: ", 
                      parseFilePaths(roots=volumes(), input$sh_start_genes)$name)
    } else {
      "Filename: "
    }
    
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
  
  sh_threads <- reactive({
    req(input$sh_threads)
    validate("Please provide the number of threads to be used when running the pipeline.")
  })
  
  mind <- reactive({
    req(input$mind)
    validate("Please provide the minimum individual to be used when running the pipeline.")
  })
  
  geno <- reactive({
    req(input$geno)
    validate("Please provide the minimum individual to be used when running the pipeline.")
  })
  
  
  observeEvent(input$sh_start,{
    # make config file
    make_short_config(in_vcf=sh_vcf_name(), out_dir=sh_out_dir(), MAF=MAF(), GQ=GQ(), 
                      DP=DP(), 
                      ibis_mt1 = ibis_mt1(), ibis_mt2=ibis_mt2())
    showNotification("Short variants prioritisation pipeline started.",
                     type="message"
    )
    Sys.sleep(3)
    system("cd scripts; ./short_variants.sh -C config/pipeline_config.config", 
           wait=FALSE)
    
  })
  #-----------------------------------------------------------------------------
  # Start SV pipeline 
  #-----------------------------------------------------------------------------
  # read in SV VCF
  shinyFileChoose(input, "sv_vcf", roots=volumes())
  sv_vcf_name <- reactive({
    if (!is.null(input$sv_vcf)){
      infile <- parseFilePaths(roots=volumes(), input$sv_vcf)$datapath
      if (length(infile) > 0){
        # check vcf has right extension (.vcf.gz)
        ifelse(substr(infile, nchar(infile)-6,
                      nchar(infile))==".vcf.gz",
               infile,
               validate("Invalid file: Input file is not a compressed VCF file (.vcf.gz).")
        )
      }
    }
  })
  
  output$sv_vcf_name <- renderText({
    req(input$sv_vcf)
    if (!is.null(input$sv_vcf)) {
      infile <- paste("Input VCF file: ", 
                      parseFilePaths(roots=volumes(), input$sv_vcf)$name)
    } else {
      "Input VCF file: "
    }
    
  })
  
  shinyDirChoose(input, "sv_outdir", roots=volumes())
  
  sv_out_dir <- reactive({
    req(input$sv_outdir)
    out_dir <- substring(parseDirPath(roots=volumes(), input$sv_outdir), 2)
  })
  
  output$sv_outdir_txt <- renderText({
    req(input$sv_outdir)
    prefix <- "Selected output folder: "
    if (!is.null(input$sv_outdir)){
      paste0(prefix, sv_out_dir()) 
    } else {
      prefix
    }
  })
  
  shinyFileChoose(input, "sv_start_genes", roots=volumes())
  
  sv_start_genes <- reactive({
    if (!is.null(input$sv_start_genes)){
      infile <- parseFilePaths(roots=volumes(), input$sv_start_genes)$datapath
      if (length(infile) > 0){
        ext=tools::file_ext(infile)
        # check segment file 
        switch(ext,
               xlsx=infile,
               txt=infile,
               validate("Invalid file: Input file is not an excel (.xlsx) or text file (.txt")
        )
      }
    }
  })
  
  output$sv_start_genes_name <- renderText({
    req(input$sv_start_genes)
    if (!is.null(input$sv_start_genes)) {
      infile <- paste("Filename: ", 
                      parseFilePaths(roots=volumes(), input$sv_start_genes)$name)
    } else {
      "Filename: "
    }
    
  })
  
  sv_threads <- reactive({
    req(input$sv_threads)
    validate("Please provide the number of threads to be used when running the pipeline.")
  })
  
  shinyFileChoose(input, "sv_start_ibis_seg", roots=volumes())
  
  sv_start_ibis_seg <- reactive({
    if (!is.null(input$sv_start_ibis_seg)){
      infile <- parseFilePaths(roots=volumes(), input$sv_start_ibis_seg)$datapath
      if (length(infile) > 0){
        ext=tools::file_ext(infile)
        # check segment file 
        switch(ext,
               seg=infile,
               validate("Invalid file: Input file is not an IBD segment file (.seg).")
        )
      }
    }
  })
  
  output$sv_ibis_seg <- renderText({
    req(input$sv_start_ibis_seg)
    if (!is.null(input$sv_start_ibis_seg)) {
      infile <- paste("IBIS IBD segment file: ", 
                      parseFilePaths(roots=volumes(), input$sv_start_ibis_seg)$name)
    } else {
      "IBIS IBD segment file: "
    }
    
  })
  
  sv_email <- reactive({
    req(input$sv_email)
    if (grep("\\w.+@\\w+.*\\..{2,4}", input$sv_email)) {
      
    } else {
      validate("Invalid email address")
    }
  })
  
  observeEvent(input$sv_start,{
    # make SV config file 
    sv_config(in_vcf=sv_vcf_name(), out_dir=sv_out_dir(), 
              sv_ibis_seg=sv_start_ibis_seg(),
              email=sv_email(), sv_threads=sv_threads(), genes=sv_start_genes())
    showNotification("Short variants prioritisation pipeline started.",
                     type="message"
    )
    Sys.sleep(3)
    system("cd scripts; ./structural_variants.sh -c config/pipeline_sv.config", 
           wait=FALSE)
    
  })
  
#-------------------------------------------------------------------------------  
# Short Variants tab
#-------------------------------------------------------------------------------
  
  # read short variants csv
  col_types <- list(CHROM='f', ID='c', REF='c', ALT='c', FILTER='f', ALLELE='f',
                    CONSEQUENCE='c', IMPACT='f', SYMBOL='f', GENE='c',
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
  shinyFileChoose(input, "short_tsv", roots=volumes())

  short_data <- reactive({
    if (!is.null(input$short_tsv)){
      infile <- parseFilePaths(roots=volumes(), input$short_tsv)$datapath
      if (length(infile) > 0){
        ext=tools::file_ext(infile)
        switch(ext,
               tsv=vroom::vroom(infile, delim="\t", 
                                col_names=TRUE, col_types = col_types),
               txt=vroom::vroom(infile, delim="\t",
                                col_names=TRUE, col_types = col_types),
               stop("Invalid file: Please upload a tsv/text file")
        )
      }
    }
   
  })
  
  output$short_tsv_name <- renderText({
    if (!is.null(input$short_tsv)) {
      infile <- paste("Variants file:", 
                      parseFilePaths(roots=volumes(), input$short_tsv)$name)
    } else {
      "Variants file:"
    }
  })
  
  output$sh_tsv_label <- renderText({
    "<b>Upload pipeline output file (.tsv/.txt)</b>"
  })
  
  shinyFileChoose(input, "ibd_seg", roots=volumes())
  ibd_data <- reactive({
    if(!is.null(input$ibd_seg)) {
      infile <- parseFilePaths(roots=volumes(), input$ibd_seg)$datapath
      if (length(infile) > 0) {
        ext=tools::file_ext(infile)
        ifelse(ext=="seg",
               infile,
               validate("Invalid file: Please upload an IBIS IBD segment file (.seg)")
        )
      }
    }
  })
  
  output$sh_seg_label <- renderText({
    "<b>Upload (IBIS) IBD segment file (.seg) from short variants pipeline</b>"
  })
  
  output$ibd_seg_name <- renderText({
    if (!is.null(input$ibd_seg)) {
      infile <- paste("IBD segment file: ", 
                      parseFilePaths(roots=volumes(), input$ibd_seg)$name)
    } else {
      "IBD segment file: "
    }
  })
  
  # render ideogram
  output$ideogram_plot <- renderIdeogram({
    if (!is.null(ibd_data())){
      ideogram({ibd_data()})
    }
    
  })
  
  # read genes list
  shinyFileChoose(input, "sh_gene_list", roots=volumes())
  sh_genes <- reactive({
    if (!is.null(input$sh_gene_list)) {
      infile <- parseFilePaths(roots=volumes(), input$sh_gene_list)$datapath
      if (length(infile) > 0) {
        ext=tools::file_ext(infile)
        switch(ext,
               xlsx=read_excel(infile, col_names=c("gene")),
               txt=vroom::vroom(infile, delim="\n",
                                col_names=c("gene")),
               stop("Invalid gene list file: Please upload a valid genes list file.")
        ) 
      }
    }
  })
  
  output$sh_gene_label <- renderText({
    "<b>Upload a list of genes of interest (.xlsx, .txt)</b>"
  })
  
  output$sh_gene_name <- renderText({
    if (!is.null(input$sh_gene_list)) {
      infile <- paste("Gene list file: ", 
                      parseFilePaths(roots=volumes(), input$sh_gene_list)$name)
    } else {
      "(Optional) Gene list file: "
    }
  })
  
  output$sh_vars_total <- renderText({
    if (!is.null(short_data())) {
      paste("Total number of variants: ", nrow(short_data()))
    } else {
      ""
    }
  })
  
  output$clinvar_vars <- renderText({
    if (!is.null(short_data())) {
      pathogenic <- short_data() %>% 
        filter(CLNSIG=="Pathogenic/Likely pathogenic") %>%
        count()
      paste('Pathogenic/ likely pathogenic (ClinVar) variants:', 
            pathogenic)
    } else {
      ""
    }
  })
  
  output$impact_vars <- renderText({
    if (!is.null(short_data())) {
      impact <- short_data() %>% 
        filter(IMPACT=="HIGH") %>%
        count()
      paste('High impact (LOF) variants:', 
            impact)
    } else {
      ""
    }
  })
  
  output$missense_vars <- renderText({
    if (!is.null(short_data())) {
      missenses <- short_data() %>% 
        filter((SIFT_CALL=="deleterious" | 
                 POLYPHEN_CALL=="probably damaging") & CADD_PHRED >= 20) %>%
        count()
      paste('Predicted functionally important missenses:', 
            missenses)
    } else {
      ""
    }
  })
  
  output$ibd_seg_total <- renderText({
    if (!is.null(ibd_data())) {
      ibd_segs <- read.table(ibd_data(), header = T)
      paste('Shared IBD segments:', 
            nrow(ibd_segs))
    } else {
      ""
    }
  })
  
  sh_gene_filter <- reactive({
    req(!is.null(input$sh_gene_check))
    if (input$sh_gene_check) {
      ibd_filter() %>% filter(filter_variables(ibd_filter()$SYMBOL, 
                                                pull(sh_genes(), gene)))
    } else {
      ibd_filter()
    }
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
                         choices=levels(cadd_bins())),
      checkboxInput("sh_gene_check", "Genes of interest", value=FALSE)
    ))
  
  # render short variants table
  output$short_tab <- renderDT({ 
    if (!is.null(short_data())) {
      DT::datatable(
        sh_gene_filter() %>%
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
          mutate(ID, ID=stringr::str_trunc(ID, width = 20)) %>%
          select(ID, RS, HGVSC, HGVSP,
                 SYMBOL, CONSEQUENCE, MAX_AF, IMPACT,
                 CADD_PHRED, POLYPHEN_CALL,
                 SIFT_CALL, CLNSIG), 
        escape=FALSE,
        rownames=FALSE) 
    }
    
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
  
  #-----------------------------------------------------------------------------
  # SV 
  #-----------------------------------------------------------------------------
  sv_col_types <- list(CHROM='f', ID='c', REF='c', ALT='c', FILTER='f', 
                        GENES='c', SVTYPE='f', SVLEN='d', END='d', 
                        CIGAR='c', CIPOS='c', CIEND='c', MATEID='c', 
                       EVENT='c', IMPRECISE='f', JUNCTION_QUAL='d')
  # sv tsv file 
  shinyFileChoose(input, "sv_tsv", roots=volumes())
  
  output$sv_tsv_label <- renderText({
    "<b>Upload SV pipeline output file (.tsv/.txt)</b>"
  })
  
  sv_data <- reactive({
    if (!is.null(input$sv_tsv)){
      infile <- parseFilePaths(roots=volumes(), input$sv_tsv)$datapath
      if (length(infile) > 0){
        ext=tools::file_ext(infile)
        switch(ext,
               tsv=vroom::vroom(infile, delim="\t", 
                                col_names=TRUE, col_types = sv_col_types),
               txt=vroom::vroom(infile, delim="\t",
                                col_names=TRUE, col_types = sv_col_types),
               stop("Invalid file: Please upload a tsv/text file")
        )
      }
    }
  })
  
  output$sv_tsv_name <- renderText({
    if (!is.null(input$sv_tsv)) {
      infile <- paste("SV file:", 
                      parseFilePaths(roots=volumes(), input$sv_tsv)$name)
    } else {
      "SV file:"
    }
  })
  
  # gene list for detecting overlaps with SV
  shinyFileChoose(input, "sv_gene_list", roots=volumes())
  sv_genes <- reactive({
    if (!is.null(input$sv_gene_list)){
      infile <- parseFilePaths(roots=volumes(), input$sv_gene_list)$datapath
      if (length(infile) > 0){
        ext=tools::file_ext(infile)
        switch(ext,
               xlsx=read_excel(infile, col_names=c("gene")),
               txt=vroom::vroom(infile, delim="\n",
                                col_names=c("gene")),
               stop("Invalid gene list file: Please upload a valid genes list file.")
        )
      }
    }
   
  })
  
  output$sv_gene_list_label <- renderText({
    "<b>(Optional) Upload a list of genes of interest (.xlsx/.txt)</b>"
  })
  
  output$sv_gene_list_name <- renderText({
    if (!is.null(input$sv_gene_list)) {
      infile <- paste("Gene list file:", 
                      parseFilePaths(roots=volumes(), input$sv_gene_list)$name)
    } else {
      "Gene list file:"
    }
  })
  
  sv_gene_filter <- reactive({
    req(!is.null(input$sv_gene_check))
    if (input$sv_gene_check) {
      sv_filters() %>% filter(GENES %in% sv_genes()$gene)
    } else {
      sv_filters()
    }
  })
  
  sv_filters <- reactive({
    if (!is.null(sv_data())){
      if (input$imprecise) {
        sv_data() %>% filter(
          filter_variables(sv_data()$CHROM, input$chrom) &
            filter_variables(sv_data()$SVTYPE, input$sv_type)
        )
      } else {
        sv_data() %>% filter(
          filter_variables(sv_data()$CHROM, input$chrom) &
            filter_variables(sv_data()$SVTYPE, input$sv_type) &
            !(IMPRECISE %in% c("TRUE"))
        )
      }
    }
  })
  
  sv_size <- reactive({
    if (!is.null(sv_data())) {
      list(min=min(sv_data()$SVLEN, 
                      na.rm=TRUE),
           max=max(sv_data()$SVLEN, 
               na.rm=TRUE))
    } else {
      list(min=0, max=0)
    }
   
  })
  
  output$sv_filters_ui <- renderUI({
    tagList(
      checkboxGroupInput("sv_type", "SV Type", 
                         choices=levels(sv_data()$SVTYPE),
                         selected=levels(sv_data()$SVTYPE)),
      checkboxGroupInput("chrom", "Chromosomes",
                         choices=levels(sv_data()$CHROM),
                         selected=levels(sv_data()$CHROM)),
      checkboxInput("imprecise", "Imprecise", value = TRUE),
      checkboxInput("sv_gene_check", "Genes of interest", value = FALSE)
    )
  })
  

  output$sv_table <- renderDT({
    DT::datatable(
      if (!is.null(sv_gene_filter())){
        sv_gene_filter()  %>%
          select(CHROM, START, END, ID, REF, ALT, SVTYPE, SVLEN, CIPOS, 
                 CIEND, CIGAR, GENES, MATEID, EVENT, IMPRECISE) %>% 
          mutate(REF, REF=stringr::str_trunc(REF, width = 20)) %>% 
          mutate(ALT, ALT=stringr::str_trunc(ALT, width = 20)) %>%
          mutate(GENES=map(GENES, process_multigenes))
      }, 
      escape=FALSE,
      rownames=FALSE,
      options=list(columnDefs = list(list(width = '10%', targets ='_all')))
      )
  })
  
  #-----------------------------------------------------------------------------
  # SV Summaries 
  #-----------------------------------------------------------------------------
  output$total_sv <- renderText({
    if (!is.null(sv_data())){
      paste0("Total number of SV calls: ", sv_data() %>% nrow())
    } else {
      ""
    }
  })
  
  output$ave_sv_len <- renderText({
    if (!is.null(sv_data())){
      paste0("Average SV length: ", pull(sv_data(), SVLEN) %>% 
               mean(na.rm=T) %>% round())
    } else {
      ""
    }
  })

  output$ave_ins_len <- renderText({
    if (!is.null(sv_data())){
      paste0("Average insertion length: ", filter(sv_data(), SVTYPE=="INS") %>% 
               pull(SVLEN) %>% mean(na.rm=T) %>% round())
    } else {
      ""
    }
  })
  
  output$ave_del_len <- renderText({
    if (!is.null(sv_data())){
      paste0("Average deletion length: ", filter(sv_data(), SVTYPE=="DEL") %>% 
               pull(SVLEN) %>% mean(na.rm=T) %>% round())
    } else {
      ""
    }
  })
  
  output$ave_dup_len <- renderText({
    if (!is.null(sv_data())){
      paste0("Average duplication length: ", 
             filter(sv_data(), SVTYPE=="DUP") %>% 
               pull(SVLEN) %>% mean(na.rm=T) %>% round())
    } else {
      ""
    }
  })
  
  output$ins_sum <- renderText({
    if (!is.null(sv_data())){
      paste0("Number of INS: ", sv_data() %>% filter(SVTYPE=="INS") %>% nrow())
    } else {
      ""
    }
  })
  
  output$del_sum <- renderText({
    if (!is.null(sv_data())){
      paste0("Number of DEL: ", sv_data() %>% filter(SVTYPE=="DEL") %>% nrow())
    } else {
      ""
    }
  })
  
  output$dup_sum <- renderText({
    if (!is.null(sv_data())){
      paste0("Number of DUP: ", sv_data() %>% filter(SVTYPE=="DUP") %>% nrow())
    } else {
      ""
    }
  })
  
  output$bnd_sum <- renderText({
    if (!is.null(sv_data())){
      paste0("Number of BND: ", sv_data() %>% filter(SVTYPE=="BND") %>% nrow())
    } else {
      ""
    }
  })
  
  output$imprecise_sum <- renderText({
    if (!is.null(sv_data())){
      paste0("Number of imprecise variants: ", sv_data() %>% 
               filter(IMPRECISE==TRUE) %>% nrow())
    } else {
      ""
    }
  })
  
  output$sv_genes_sum <- renderText({
    if (!is.null(sv_data())){
      paste0("Number of SVs overlapping genes: ", sv_data() %>% 
               filter(GENES!=".") %>% nrow())
    } else {
      ""
    }
  })
  
  output$max_sv_len <- renderText({
    if (!is.null(sv_data())){
      paste0("Maximum SV length: ", sv_data() %>% 
               pull(SVLEN) %>% max(na.rm=T))
    } else {
      ""
    }
  })
  
  output$genes_list_sum <- renderText({
    if (!is.null(sv_genes())){
      paste0("Number of SV overlapping with genes of interest: ", sv_data() %>% 
               filter(GENES %in% sv_genes()$gene) %>% nrow())
    } else {
      ""
    }
  })
  
  # download variants results handler
  output$sv_download <- downloadHandler(
    filename=function(){
      paste0(
        tools::file_path_sans_ext(input$sv_tsv$name), "_filtered.tsv")
    },
    content=function(file){
      write.table({sv_filters()}, file)
    })
  
}
