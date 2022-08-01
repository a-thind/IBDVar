# Shiny App
#

library(shiny)
library(shinydashboard)
library(shinyFiles)
library(biomaRt)
library(ideogram)
library(DT)

# function to create select inputs for factor variables
filters_ui <- function(var) {
  levels = levels(var)
  selectInput(var, var, choices=levels, selected = NULL, multiple = TRUE)
}

# make_link <- function(gene_id) {
#   sprintf(
#     "<a href=https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s>%s</a>",
#     gene_id, gene_id)
# }

# source, output folder and config
#

# Define UI for application
# sidebar
sidebar <- dashboardSidebar(
  sidebarMenu(
    id="sidebar_id",
    menuItem("Start Pipeline", tabName = "pipeline"),
    menuItem("Short Variants", tabName="short_vars"),
    menuItem("Structural Variants", tabName="sv_tab"),
# Filters side panel    
#-------------------------------------------------------------------------------
    conditionalPanel(
      'input.sidebar_id == "short_vars"',
      tags$h3("Filters"),
      uiOutput("Filters"),
      sliderInput("cadd_filter", "CADD Score", value=20, min=1, max=100),
      checkboxGroupInput("impact_filter", "VEP Impact",
                         choices=c("high", "noderate", "low")),
      checkboxGroupInput("sift_filter", "SIFT",
                         choices=c("deleterious", 
                                   "deleterious low confidence",
                                   "tolerated",
                                   "tolerated low confidence")),
      checkboxGroupInput("polyphen_filter", "PolyPhen", 
                         choices=c("probably damaging",
                                   "possibly damaging",
                                   "benign"))
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
            column(8,
              fileInput("sh_vcf", "Upload input short variants VCF file (compressed)",
                             accept=c("vcf.gz"), width="75%"),
              shinyDirButton("sh_outdir", "Select output folder",
                                    title="Output folder", multiple=F),
              textInput("email", "Enter an email address:", width="75%"),
              actionButton("sh_start", "Start")
              )
            ),

            box(title="Advanced Configuration", collapsible = TRUE,
                collapsed = TRUE,
                status = "warning")

      ),
      fluidRow(
        box(
          title="Structural Variants",
          status = "primary",
          fileInput("sv_vcf", "Upload input structural variants VCF file (compressed).",
                    accept=c("vcf.gz"), width="50%")
        )
      )

    ),
# Short variants tab
#-------------------------------------------------------------------------------
    tabItem(
      tabName="short_vars",
      fluidRow(
        box(
          title = "IBD regions Ideogram",
          status="primary",
          width = 8,
          height="600px",
          box(
            width = 12,
            height="535px",
            ideogramOutput("ideogram_plot")
          )

        ),
        box(
          title="Files",
          width=4,
          fileInput("short_tsv", "Upload pipeline output file (.tsv/.txt)",
                    accept=c(".tsv", ".txt")),
          fileInput("ibd_seg", "Upload IBIS IBD segment file (.seg)",
                    accept=c(".seg")),
          status = "primary"
        )
      ),
      fluidRow(
        box(
          width=12,
          height="100%",
          actionButton("reset_short_tab", "Reset"),
          downloadButton("download", "Download"),
          tags$br(),
          DTOutput("short_tab"),
          status = "primary"
        )
      )
    ), 
    # SV tab
    #-----------------
    tabItem(
      tabName="sv_tab", 
      fluidRow(
        box(
          title = "IBD regions Ideogram",
          status="primary",
          width = 8,
          height="600px",
          box(
            width = 12,
            height="535px",
            ideogramOutput("sv_ideogram_plot")
          )
          
        ),
        box(
          title="Files",
          width=4,
          fileInput("sv_tsv", "Upload SV pipeline output file (.tsv/.txt)",
                    accept=c(".tsv", ".txt")),
          fileInput("sv_ibd_seg", "Upload IBIS IBD segment file (.seg)",
                    accept=c(".seg")),
          status = "primary",
      ))
  )
)
)

# define user interface
ui <- dashboardPage(
  dashboardHeader(title="IBD Variants"),
  sidebar,
  body
)
