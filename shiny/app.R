#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

setwd("~/Repos/MSPtoDB")
#source("lib/blobfunctions.R")
source("lib/msp2db.R")


# Define UI for application that draws a histogram
options(shiny.maxRequestSize = 30*1024^2)

ui <- fluidPage(
    #MSPtoDB(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, DBoutput, Source)
    
    # Application title
    titlePanel("Library Maker"),
    fileInput("LibInput", "Input File"),
    selectInput("FragInput", "Fragmentation mode", c("HCD", "CID")),
    numericInput("CeInput", "Collision Energy", 35, min =0, max = 00),
    selectInput("MassAnalyzerInput", "Mass Analyzer", c("IT", "OT")),
    selectInput("TMTProInput", "Add TMTPro?", c("TRUE", "FALSE")),
    selectInput("SourceInput", "Source", c("Prosit", "SpectraST")),
    textInput("DbInput", ".db file name"),
    actionButton("button", "Go"),
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
    DBoutput = reactive(input$DbInput)
    
    observeEvent(input$button, {
        outputDb<-(MSPtoDB(Library()$datapath, FragmentationMode(), MassAnalyzer(), as.character(CollisionEnergy()), TMTPro(), DBoutput(), Source()))
        outputDb
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
