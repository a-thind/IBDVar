# Shiny App GUI code

# Define UI for application
# sidebar
sidebar <- dashboardSidebar(
  sidebarMenu(
    id="sidebar_id",
    menuItem("Start Pipeline", tabName = "pipeline"),
    menuItem("Short Variants", tabName="short_vars"),
    menuItem("Structural Variants", tabName="sv_tab"),
    # Filters side panel    
    #---------------------------------------------------------------------------
    conditionalPanel(
      'input.sidebar_id == "short_vars"',
      uiOutput("filters")
    ),
    conditionalPanel(
      'input.sidebar_id == "sv_tab"',
      uiOutput("sv_filters_ui")
    )
  )
)
# dashboard body
body <- dashboardBody(
  tabItems(
    # Starting pipeline
    #-------------------------------------------------------------------------------
    tabItem(
      tabName="pipeline",
      fluidRow(
        box(title="Short Variants",
            status="primary",
            width=9,
            height = "100%",
            tags$head(
              tags$style(
                HTML(".form-control.shiny-bound-input{
                     width: 75px;
                }
                #sh_email{
                  width: 75%
                }
                #sv_email{
                  width: 75%
                }
                     "))),
              column(6,
                    tags$h4(
                      "General Settings"
                    ),
                    fluidRow(
                             column(12,
                                    htmlOutput("sh_vcf_label")
                                    )
                             ),
                    fluidRow(
                      column(3,  
                             shinyFilesButton("sh_vcf", 
                                              "Browse...",
                                              title="Upload input short variants VCF file (compressed)",
                                              filetype=c("vcf.gz"), multiple=F),
                             
                             ),
                      column(9,
                             textOutput("sh_vcf_name"),
                      )),
                      fluidRow(
                        column(12, 
                               br(),
                               htmlOutput("sh_outdir_label")),
                      ),
                      fluidRow(
                        column(3,
                          shinyDirButton("sh_outdir", "Browse...",
                          title="Output folder", multiple=F)),
                        column(9,
                          textOutput("sh_outdir_txt")
                        )
                      ),       
                    numericInput("sh_threads", "Number of threads (CPU)", value=4,
                                 min=1, max=99),
                    textInput("sh_email","Email address"),
                    br(),
                    tags$h4(
                      "PLINK Dataset Generation"
                    ),  
                    numericInput("min_AF", "Min. Allele Frequency threshold for filtering out variants for PLINK dataset", value=0.05, max=1, 
                                 min=0.01, step=0.01),
                    tags$h4("Missing genotype rates (for PLINK dataset)"),
                    numericInput("mind", "Exclude SAMPLES with more than the provided value percentage of missing call rates e.g. 0.1 excludes samples with missing call rates > 10%", value=0.05, max=1, 
                                 min=0.1, step=0.01),
                    numericInput("geno", "Exclude VARIANTS with missing call rates above the provided value: e.g. 0.1 excludes variants with missing call rates > 10%", value=0.01, max=1, 
                                 min=0.1, step=0.01)
              ),
              column(6,
                     tags$h4(
                       "QC Filtering"
                     ),
                     numericInput("GQ", "Minimum genotype quality (GQ) per sample threshold", value=20,
                                  min=1, max=99),
                     numericInput("DP", "Minimum read depth (DP) per sample threshold", value=10, min=1,
                                  max=100),
                     br(),
                     tags$h4(
                       "IBD Segment Detection"
                     ),
                     numericInput("ibis_mt1", 
                                  "Minimum number of (SNP) markers to call a region IBD1",
                                  value=50, min=10, max=1000),
                     numericInput("ibis_mt2",
                                  "Minimum number of (SNP) markers to call region IBD2",
                                  value=10, min=1, max=400),
                     br(),
                     tags$h4("Protein-impacting Variant Selection"),
                     numericInput("MAF", "Max Allele Frequency in any of the following populations: gnomAD, 1000 genomes or ESP", value=0.05, max=1, 
                                  min=0.01, step=0.01),
                     actionButton("sh_start", "Start")
                     )
            
        ),
        box(
          title="Structural Variants",
          status = "primary",
          height="100%",
          width=3,
          fluidRow(
            column(12,
                   htmlOutput("sv_vcf_label")
                   )
          ),
          fluidRow(
            column(3,
               shinyFilesButton("sv_vcf", "Browse...", 
                  title="Upload input structural variants VCF file (compressed)",
                  filetype=c(".vcf.gz"), multiple=F)  
            ),
            column(6,
                   
            )
          ),
          fluidRow(
            column(3,
                   shinyDirButton("sv_outdir", "Select output folder",
                                  title="Output folder", multiple=F)
                   ),
            column(6,
                   textOutput("sv_outdir_txt")
                   )
            ),
          fluidRow(
            column(3,
                   shinyFilesButton("sv_start_ibis_seg", "Browse...", 
                                    title="IBIS IBD segment file (.seg)",
                                    filetype=c(".seg"), multiple=F)
                   ),
            column(6,
                   textOutput("sv_ibis_seg")
                   )
          ),
          column(12,
                 numericInput("sv_threads", "Number of threads", value=4,
                              min=1, max=99),
                 actionButton("sv_start", "Start")
          )
        )
      )
    ),
    #---------------------------------------------------------------------------
    # Short variants tab
    #---------------------------------------------------------------------------
    tabItem(
      tabName="short_vars",
      fluidRow(
        box(
          title = "IBD regions Ideogram",
          status="primary",
          width = 9,
          height="610px",
          box(
            width = 12,
            height="540px",
            ideogramOutput("ideogram_plot")
          )
        ),
        column(
          width=3,
          box(
            title="Files",
            width=12,
            fluidRow(column(12,
                            htmlOutput("sh_tsv_label")
            )),
            fluidRow(
              column(3,
                     br(),
                     shinyFilesButton(
                       "short_tsv", "Browse...", 
                        title="Upload pipeline output file (.tsv/.txt)",
                        filetype=c(".tsv", ".txt"), multiple = FALSE)),
              column(9,
                     br(),
                     textOutput("short_tsv_name"))
            ),
            fluidRow(
              column(12,
                     br(),
                     htmlOutput("sh_seg_label")
              )
            ),
            fluidRow(
              column(3,
                     br(),
                     shinyFilesButton(
                       "ibd_seg", "Browse...", 
                        title="Upload (IBIS) IBD segment file (.seg)",
                        filetype=c(".tsv", ".txt"), 
                        multiple = FALSE)
              ),
              column(9,
                     br(),
                     textOutput("ibd_seg_name"))
            ),
            fluidRow(
              column(12,
                     br(),
                     htmlOutput("sh_gene_label")
              )),
            fluidRow( 
              column(3,
                     br(),
                     shinyFilesButton("sh_gene_list", "Browse...", 
                                      title="(Optional) Upload a list of genes of interest (.xlsx, .txt)",
                                      filetype=c(".xlsx", ".txt"), 
                                      multiple=FALSE)
                     ),
              column(9,
                     br(),
                     textOutput("sh_gene_name"))
            ),
            status = "primary"
          ),
          box(
            title="Summary",
            width=12,
            height="250px",
            status="primary",
            column(12,
              fluidRow(
                textOutput("sh_vars_total"),
                br()
              ),
              fluidRow(
                textOutput("clinvar_vars"),
                br()
              ),
              fluidRow(
                textOutput("impact_vars"),
                br()
              ),
              fluidRow(
                textOutput("missense_vars"),
                br()
              ),
              fluidRow(
                textOutput("ibd_seg_total")
              )
            ))
        )
      ),
      fluidRow(
        box(
          width=12,
          height="100%",
          column(12,
                 fluidRow(
                   downloadButton("download", "Download")
                 ),
                 br(),
                 fluidRow(
                   DTOutput("short_tab")
                 )
                )
          ,
          status = "primary"
        )
      )
    ), 
    #---------------------------------------------------------------------------
    # SV tab
    #---------------------------------------------------------------------------
    tabItem(
      tabName="sv_tab", 
      fluidRow(
        box(
          title = "Summary",
          status="primary",
          width = 9,
          height="100%",
          box(
            width = 12,
            height="100%"
          )
        ),
        box(
          title="Files",
          width=3,
          fluidRow(
            column(12,
                   htmlOutput("sv_tsv_label"),
                   br()
            ) 
          ), 
          fluidRow(
            column(3,
                  shinyFilesButton("sv_tsv", "Browse...",title="Upload SV pipeline output file (.tsv/.txt)",
                                   filetype=c(".tsv", ".txt"), multiple=F)
                  ),
            column(9,
                   textOutput("sv_tsv_name")
              )
          ),
          fluidRow(
            column(12,
                   htmlOutput("sv_gene_list_label"),
                   br()
            )
          ),
          fluidRow(
            column(3,
                   shinyFilesButton("sv_gene_list", "Browse...",
                                    title="Upload a list of genes of interest (.xlsx/.txt)", 
                                    filetype=c(".xlsx", ".txt"), multiple=F)
            ),
            column(9,
                   textOutput("sv_gene_list_name")
            )
          ),
          status = "primary",
          textOutput("text")
        )),
      fluidRow(
        box(
          title = "Stuctural Variants",
          status="primary",
          height="100%",
          width = "100%",
          downloadButton("sv_download", "Download"),
          tags$br(),
          DTOutput("sv_table")
        )
      )
    )
  )
)

# define user interface
ui <- dashboardPage(
  dashboardHeader(title="IBD Variants"),
  sidebar,
  body
)
