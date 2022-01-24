
library(shiny)
library(shinybusy)
library(shinydashboard)
library(RSQLite)                                                                                                                                                                                          
library(stringr)                                                                                                                                                                                          
library(vroom)                                                                                                                                                                                             
library(stringi)                                                                                                                                                                                          
library(tidyverse)                                                                                                                                                                                            
library(readr)                                                                                                                                                                                            
library(data.table)                                                                                                                                                                                       
library(purrr)                                                                                                                                                                                            
library(blob)
library(BiocParallel)
library(shinyjs)

setwd("~/Repos/MSPtoDB")
source("lib/msp2db.R")
source("lib/LibraryParserv2.R")


# Define UI for application that draws a histogram
options(shiny.maxRequestSize = 3000*1024^2)

#set port and host options
options(shiny.port = 3838)
options(shiny.host = "0.0.0.0")

header <- dashboardHeader(
  title = "Library Maker"
)

side <- dashboardSidebar(
  sidebarMenu(
    menuItem("Library Maker", tabName = "maker", icon = icon("font-awesome"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "maker",
            fluidPage(
              useShinyjs(),
              titlePanel("Spectral Library to mzVault .db Converter"),
              fileInput("LibInput", "Input Files", accept=c('.msp','.MSP','.sptxt','.blib'), multiple = TRUE),
              selectInput("FragInput", "Fragmentation", c("Read From File"= '', "HCD", "CID")),
              textInput("CeInput", "Normalized Collision Energy (NCE)","", placeholder =  "Read from file" ),
              selectInput("MassAnalyzerInput", "Mass Analyzer", c("OT", "IT")),
              #selectInput("TMTProInput", "Add TMTPro?", c("FALSE", "TRUE")),
              #selectInput("Filter", "Filter peaks?", c("Yes", "No")),
              numericInput("massoffset", "Mass Offset: ", 0),
              checkboxInput("Filter", "Filter peaks?", value = FALSE, width = NULL),
              hidden(checkboxGroupInput(
                "IonTypes",
                "Ion types included?",
                choices = c("a","b","c","y","z"),
                selected = NULL,
                inline = FALSE,
                width = NULL,
                choiceNames = NULL,
                choiceValues = NULL
              )),
              textOutput("text"),
              hidden(numericInput("topX", "Top N Peaks only", 150)),
              hidden(numericInput("cutoff", "% intensity cutoff", 0)),
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
  TMTPro = reactive(input$TMTProInput)
  DBoutput = reactive(input$DbInput)
  Filter = reactive(input$Filter)
  topX = reactive(input$topX)
  cutoff = reactive(input$cutoff)
  observe({
    toggle(id="topX", condition = input$Filter)})
  observe({
    toggle(id="cutoff", condition = input$Filter)})
  observe({
      toggle(id="IonTypes", condition = input$Filter)})
  output$text <- renderText({ paste(input$IonTypes)})
  output$downloadData <- downloadHandler(
   # filename = function() { paste('library-',format(Sys.time(), "%Y-%m-%d_%I-%p"),'.db',sep='') },
    filename = function() { paste0(substring(Library()$name,0,nchar(Library()$name)-3), 'db') },
   # content = function(outFile) { DBbuilder(Library()$datapath,FragmentationMode(), MassAnalyzer(), as.character(CollisionEnergy()), TMTPro(), "test.db",FALSE,Source(), topX(), cutoff()) }
    content= function(x) {
      DBbuilder(Library()$datapath, "CID", "OT", "35", FALSE ,FALSE , x, 150, 0)
    } )
    
}

# Run the application 
shinyApp(ui = ui, server = server)