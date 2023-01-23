
library(shiny)
library(shinybusy)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinyWidgets)

setwd("~/Repos/MSPtoDB")
source("lib/msp2db.R")
source("lib/LibraryParserv2.R")


# Define UI for application that converts spectral library formats
options(shiny.maxRequestSize = 10000*1024^2)

#set port and host options
options(shiny.port = 3838)
options(shiny.host = "0.0.0.0")
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5000)

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
              textOutput("var"),            
              fileInput("LibInput", "Input Files", accept=c('.msp','.MSP','.sptxt','.blib'), multiple = TRUE),
              selectInput("FragInput", "Fragmentation", c("Read From file"= '', "HCD", "CID")),
              textInput("CeInput", "Normalized Collision Energy (NCE)","Read from file", placeholder =  "Read from file" ),
              selectInput("MassAnalyzerInput", "Mass Analyzer", c("FT", "IT")),
              fileInput("massOffset", "Mass Offset CSV File: ", accept=c(".csv")),
              fileInput("CompClassInput", "Compound Class CSV File: ", accept=c(".csv")),
              div(switchInput("Filter", label="Filter", value = FALSE),
                  style = "font-size: 20px !important; text-align:left;"
              ),
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
              textOutput("txt"),
              hidden(switchInput("TMTInput",label= "Remove TMT-ions?",value= FALSE)),
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


# Define server logic required to create a spectral library
server <- function(input, output) {
  output$var <- renderText({ 
    grepl("blib",Library()$name[1])
  })
  Library = reactive(input$LibInput)
  FragmentationMode = reactive(input$FragInput)
  MassAnalyzer =  reactive(input$MassAnalyzerInput)
  CollisionEnergy = reactive(input$CeInput)
  TMTPro = reactive(input$TMTInput)
  Filter = reactive(input$Filter)
  Decoy = reactive(input$Decoy)
  DBoutput = reactive(input$DbInput)
  topX = reactive(input$topX)
  cutoff = reactive(input$cutoff)  
  rtalign = reactive(input$RTinput)
  output$text <- renderText({ Library()$name }) 
  IonTypes=reactive(input$IonTypes)
  observe({
    toggle(id="topX", condition = input$Filter)})
  observe({
    toggle(id="cutoff", condition = input$Filter)})
  observe({
      toggle(id="IonTypes", condition = input$Filter)})
  output$downloadData <- downloadHandler(
   # filename = function() { paste('library-',format(Sys.time(), "%Y-%m-%d_%I-%p"),'.db',sep='') },
    filename = function() { paste0(gsub("[^.]+$", "", Library()$name), 'db') },
     content= function(x) {
       Y<-isolate({
         massOff<-input$massOffset
         top<-input$topX
         cutoff<-input$cutoff
         IonTypes<-input$IonTypes
         Filter<-input$Filter
         FragmentationMode = (input$FragInput)
         MassAnalyzer =  (input$MassAnalyzerInput)
         CollisionEnergy = (input$CeInput)
         CompoundClassArg = input$CompClassInput
         TMTPro = (input$TMTInput)
         })
        DBbuilder(Library=Library(), FragmentationMode=FragmentationMode, MassAnalyzer=MassAnalyzer(), CollisionEnergy=CollisionEnergy, CompoundClass=CompoundClassArg,
                           Filter=Filter, DBoutput=x, topX=top, TMTPro = TMTPro, cutoff=cutoff, massOffset=massOff, IonTypes=IonTypes)
   
    } )
  output$txt <-  renderText({
    input$IonTypes
  })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
