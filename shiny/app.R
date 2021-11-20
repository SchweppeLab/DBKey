#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinybusy)
library(shinydashboard)
library(RSQLite)                                                                                                                                                                                          
library(stringr)                                                                                                                                                                                          
library(Rcpp)                                                                                                                                                                                             
library(stringi)                                                                                                                                                                                          
library(dplyr)                                                                                                                                                                                            
library(readr)                                                                                                                                                                                            
library(data.table)                                                                                                                                                                                       
library(purrr)                                                                                                                                                                                            
library(blob)

setwd("~/Repos/MSPtoDB")
source("lib/blobfunctions.R")
source("lib/msp2db.R")


# Define UI for application that draws a histogram
options(shiny.maxRequestSize = 300*1024^2)

#set port and host options
options(shiny.port = 3838)
options(shiny.host = "0.0.0.0")

header <- dashboardHeader(
  title = "Library Maker"
)

side <- dashboardSidebar(
          sidebarMenu(
             menuItem("Library Maker", tabName = "maker", icon = icon("dashboard"))
          )
      )

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "maker",
      fluidPage(
        titlePanel("Spectral Library to mzVault .db Converter"),
		fileInput("LibInput", "Input Files", accept=c('.msp','.MSP','.sptxt'), multiple = TRUE),
        selectInput("FragInput", "Fragmentation mode", c("HCD", "CID")),
        numericInput("CeInput", "Normalized Collision Energy (NCE)", 35, min =0, max = 00),
        selectInput("MassAnalyzerInput", "Mass Analyzer", c("OT", "IT")),
        selectInput("TMTProInput", "Add TMTPro?", c("FALSE", "TRUE")),
        selectInput("SourceInput", "Source", c("Prosit", "SpectraST")),
        numericInput("topX", "Top N Peaks only", 150),
        numericInput("cutoff", "% intensity cutoff", 0),
        downloadButton(label = "Generate and Download .db", "downloadData"),
        add_busy_spinner(spin = "fading-circle", margins = c(80, 300), position='top-left',height="200px",width="200px")
      )
    )
  )
 )

ui <- dashboardPage(
        header,
        side,
        body
        )


# Define server logic required to draw a histogram
server <- function(input, output) {
    Library = reactive(input$LibInput)
    FragmentationMode = reactive(input$FragInput)
    MassAnalyzer =  reactive(input$MassAnalyzerInput)
    CollisionEnergy = reactive(input$CeInput)
    Source = reactive(input$SourceInput)
    TMTPro = reactive(input$TMTProInput)
    DBoutput = reactive(input$DbInput)
    topX = reactive(input$topX)
    cutoff = reactive(input$cutoff)
    output$downloadData <- downloadHandler(
        filename = function() { paste('library-',format(Sys.time(), "%Y-%m-%d_%I-%p"),'.db',sep='') },
        content = function(outFile) { MSPtoDB(Library(),FragmentationMode(), MassAnalyzer(), as.character(CollisionEnergy()), TMTPro(), outFile, Source(), topX(), cutoff()) }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
