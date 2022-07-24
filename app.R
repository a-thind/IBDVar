#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(biomaRt)
library(ideogram)

csqs <- c("missense", "frame shift", "stop gain", "stop loss")
vep_impact <- c("high", "moderate", "low")

make_link <- function(gene_id) {
  sprintf(
    "<a href=https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s>%s</a>", 
    gene_id, gene_id)
} 

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("IBD Variant Analysis"),
    tabsetPanel(
      id="tabset",
      tabPanel("Start Pipeline", 
               h2("Import Data"),
               fileInput("short_vcf", "Short variants VCF file", 
                         buttonLabel="Upload",
                         accept=c("vcf", "vcf.gz")),
               fileInput("short_md5sum", "Md5Sum file", buttonLabel="Upload"),
               tabPanel("Configuration")
      ),
      # short variants tab panel
      tabPanel( "Short Variants",
                sidebarLayout(
                  sidebarPanel(
                    tabPanel("short_files", h4("Import Data"), 
                             fileInput("short_csv", 
                                       "Short variants output file (csv)", 
                                       buttonLabel="Upload"),
                             fileInput("ibd_seg", "IBD segment (IBIS) file", 
                                       buttonLabel="Upload",
                                       accept=c("seg"))
                    )
                    ,
                    hr(),
                    checkboxGroupInput("csqs", h4("Consequence"), csqs),
                    hr(),
                    checkboxGroupInput("vep_impact",h4("VEP Impact"), vep_impact)
                  ),
                  
                  mainPanel(
                    fluidRow(ideogramOutput("ideogram"), 
                                    dataTableOutput("short_vars"))
                  )
                )
                
      ), 
      tabPanel("SV",
               sidebarLayout(
                 sidebarPanel(
                   fileInput("short_csv",
                             "SV pipeline output file (csv)", 
                             buttonLabel="Upload"),
                   fileInput("ibd_seg", "IBD segment (IBIS) file", 
                             buttonLabel="Upload",
                             accept=c("seg"))
                 ), 
                 mainPanel()
               ))
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # reading short csv data
  short_data <- reactive({
    req(input$short_csv)
    # check file exts
    ext=tools::file_ext(input$short_csv$name)
    switch(
      ext,
      csv=vroom::vroom(input$short_csv$datapath, delim=","),
      tsv=vroom::vroom(input$short_csv$datapath, delim="\t"),
      txt=vroom::vroom(input$short_csv$datapath, delim="\t"),
      validate("Invalid short variants file, please upload a csv or tsv")
    )
    
  })
  
  ibd_data <- reactive({
    req(input$ibd_seg)
    ext=tools::file_ext(input$ibd_seg$name)
    switch(
      ext,
      seg=input$ibd_seg$datapath,
      validate("Invalid IBD segment file, please upload IBIS IBD segment file")
    )
  })
  
  output$ideogram <- renderIdeogram({
    ideogram(ibd_data())
    })
  
  # render short variant table
  output$short_vars <- renderDataTable(
    {short_data()}, options=list(pageLength=10), escape = FALSE)
}

# Run the application 
shinyApp(ui = ui, server = server)
