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
#source("lib/blobfunctions.R")
source("lib/msp2db.R")


# Define UI for application that draws a histogram
options(shiny.maxRequestSize = 300*1024^2)

#set port and host options
options(shiny.port = 3838)
options(shiny.host = "0.0.0.0")

ui <- fluidPage(
    #MSPtoDB(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, DBoutput, Source, topX, cutoff)
    
    # Application title
    titlePanel("Library Maker"),
    fileInput("LibInput", "Input File"),
    selectInput("FragInput", "Fragmentation mode", c("HCD", "CID")),
    numericInput("CeInput", "Collision Energy", 35, min =0, max = 00),
    selectInput("MassAnalyzerInput", "Mass Analyzer", c("IT", "OT")),
    selectInput("TMTProInput", "Add TMTPro?", c("TRUE", "FALSE")),
    selectInput("SourceInput", "Source", c("Prosit", "SpectraST")),
    numericInput("topX", "Top N Peaks only", 150),
    numericInput("cutoff", "% intensity cutoff", 0),
    downloadButton(label = "Generate and Download .db", "downloadData"),
    add_busy_spinner(spin = "fading-circle")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    Library = reactive(input$LibInput)
    FragmentationMode = reactive(input$FragInput)
    MassAnalyzer =  reactive(input$MassAnalyzerInput)
    CollisionEnergy = reactive(input$CeInput)
    Source = reactive(input$SourceInput)
    TMTPro = reactive(input$TMTProInput)
    topX = reactive(input$topX)
    cutoff = reactive(input$cutoff)
    output$downloadData <- downloadHandler(
        filename = function() { paste('library-',Sys.Date(),'.db',sep='') },
        content = function(outFile) { MSPtoDB(Library()$datapath,FragmentationMode(), MassAnalyzer(), as.character(CollisionEnergy()), TMTPro(), outFile, Source(), topX(), cutoff()) }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
