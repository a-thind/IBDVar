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
            column(8,
                   fileInput("sh_vcf", "Upload input short variants VCF file (compressed)",
                             accept=c("vcf.gz"), width="75%"),
                   shinyDirButton("sh_outdir", "Select output folder",
                                  title="Output folder", multiple=F),
                   textOutput("sh_outdir_txt"),
                   numericInput("GQ", "Minimum genotype quality (GQ) per sample threshold", value=20,
                                min=1, max=99),
                   numericInput("DP", "Minimum read depth (DP) per sample threshold", value=10, min=1,
                                max=100),
                   numericInput("MAF", "Minor Allele Frequency (MAF) in any of the following populations: gnomAD, 1000 genomes or ESP", value=0.05, max=1, 
                                min=0.01, step=0.01),
                   numericInput("ibis_mt1", 
                                "Minimum number of (SNP) markers to call a region IBD1",
                                value=50, min=10, max=1000),
                   numericInput("ibis_mt2",
                                "Minimum number of (SNP) markers to call region IBD2",
                                value=10, min=1, max=400),
                   numericInput("mind", "Exclude samples with more than the provided value percentage of missing genotype data e.g. 0.1 excludes samples with > 10% missing genotype data", value=0.05, max=1, 
                                min=0.1, step=0.01),
                   numericInput("geno", "Selects variants with missing call rates lower than the provided value to be removed", value=0.01, max=1, 
                                min=0.1, step=0.01),
                   numericInput("sh_threads", "Number of threads", value=4,
                                min=1, max=99),
                   actionButton("sh_start", "Start")
            )
        ),
        box(
          title="Structural Variants",
          status = "primary",
          height="100%",
          column(8,
                 fileInput("sv_vcf", "Upload input structural variants VCF file (compressed)",
                           accept=c(".gz"), width="75%"),
                 shinyDirButton("sv_outdir", "Select output folder",
                                title="Output folder", multiple=F),
                 fileInput("sv_start_ibis_seg", 
                           "IBIS IBD segment file (.seg)",
                           accept=c(".seg"), width="75%"),
                 textOutput("sv_outdir_txt"),
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
          height="600px",
          box(
            width = 12,
            height="535px",
            ideogramOutput("ideogram_plot")
          )
        ),
        box(
          title="Files",
          width=3,
          fluidRow(column(12,
                          htmlOutput("sh_tsv_label")
          )),
          fluidRow(
            column(3,
                   shinyFilesButton("short_tsv", "Browse...", 
                                    title="Upload pipeline output file (.tsv/.txt)",
                                    filetype=c(".tsv", ".txt"), multiple = FALSE)),
            column(9,
                   textOutput("short_tsv_name"))
          ),
          fluidRow(
            column(12,
                   htmlOutput("sh_seg_label")
            )
          )
          ,
          fluidRow(
            column(3,
                   shinyFilesButton("ibd_seg", "Browse...", 
                                    title="Upload (IBIS) IBD segment file (.seg)",
                                    filetype=c(".tsv", ".txt"), multiple = FALSE)
            ),
            column(9,
                   textOutput("ibd_seg_name"))
          ),
          fluidRow(
            column(12,
                    htmlOutput("sh_gene_label")
          )),
          fluidRow( 
            column(3,
                    shinyFilesButton("sh_gene_list", "Browse...", 
                    title="Upload a list of genes of interest (.xlsx, .txt)",
                    filetype=c(".xlsx", ".txt"), multiple=FALSE)),
             column(9,
                    textOutput("sh_gene_name"))
          ),
          status = "primary"
        )
      ),
      fluidRow(
        box(
          width=12,
          height="100%",
          downloadButton("download", "Download"),
          actionButton("show_all", "Show All"),
          tags$br(),
          DTOutput("short_tab"),
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
          width = 8,
          height="100%",
          box(
            width = 12,
            height="100%"
          )
        ),
        box(
          title="Files",
          width=4,
          fileInput("sv_tsv", "Upload SV pipeline output file (.tsv/.txt)",
                    accept=c(".tsv", ".txt"), width="100%"),
          fileInput("sv_gene_list", 
                    "Upload a list of genes of interest (.xlsx/.txt)", 
                    width="75%", accept=c(".xlsx", ".txt")),
          status = "primary",
          textOutput("text")
        )),
      fluidRow(
        box(
          title = "SV Overlap",
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
