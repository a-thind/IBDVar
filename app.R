#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

csqs <- c("missense", "frame shift", "stop gain", "stop loss")
vep_impact <- c("high", "moderate", "low")

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
                    dataTableOutput("short_vars")
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
  short_data <- reactive({
    req(input$short_csv)
    read.table(input$short_csv$datapath)
  })
  output$short_vars <- renderDataTable(
    {short_data()}, options=list(pageLength=5))
}

# Run the application 
shinyApp(ui = ui, server = server)
