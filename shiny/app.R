
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
              actionButton("previewButton", "Preview Entry"),
              #plotOutput("previewPlot"),
              hidden(numericInput("topX", "Top N Peaks only", 150)),
              hidden(numericInput("cutoff", "% intensity cutoff", 0)),
              div(switchInput("AdjustFragments", label="Shift Fragments", value = FALSE),
                  style = "font-size: 20px !important; text-align:left;"
              ),
              hidden(numericInput(
                "oldMod",
                "Current modification mass",0
              )),
              hidden(numericInput(
                "newMod",
                "New modification mass",0
              )),
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
  db_created <- reactiveVal(NULL)
  Library = reactive(input$LibInput)
  FragmentationMode = reactive(input$FragInput)
  MassAnalyzer =  reactive(input$MassAnalyzerInput)
  CollisionEnergy = reactive(input$CeInput)
  Filter = reactive(input$Filter)
  Decoy = reactive(input$Decoy)
  DBoutput = reactive(input$DbInput)
  topX = reactive(input$topX)
  cutoff = reactive(input$cutoff)  
  rtalign = reactive(input$RTinput)
  output$text <- renderText({ Library()$name }) 
  IonTypes=reactive(input$IonTypes)
  adjustfragments = reactive(input$AdjustFragments)
  oldmod = reactive(input$oldMod)
  newmod = reactive(input$newMod)
  observe({
    toggle(id="topX", condition = input$Filter)})
  observe({
    toggle(id="cutoff", condition = input$Filter)})
  observe({
      toggle(id="IonTypes", condition = input$Filter)})
  observe({
    toggle(id="oldMod", condition = input$AdjustFragments)})
  observe({
    toggle(id="newMod", condition = input$AdjustFragments)})
  # observeEvent(input$downloadButton, {
  #   # get a list of all input names
  #   input_names <- names(input)
  #   
  #   # create a list of all input values
  #   input_values <- lapply(input_names, function(x) {
  #     input[[x]]
  #   })
  #   
  #   # print the list of input values
  # })
  # get_random_entry <- function() {
  #   # Check if the input is available
  #   if (is.null(Library())) {
  #     return(NULL)
  #   }
  #   if (is.null(db_created())) {
  #   # Get the reactive values
  #   FragmentationMode <- isolate(input$FragInput)
  #   MassAnalyzer <- isolate(input$MassAnalyzerInput)
  #   CollisionEnergy <- isolate(input$CeInput)
  #   CompoundClassArg <- isolate(input$CompClassInput)
  #   Filter <- isolate(input$Filter)
  #   top <- isolate(input$topX)
  #   cutoff <- isolate(input$cutoff)
  #   massOff <- isolate(input$massOffset)
  #   IonTypes <- isolate(input$IonTypes)
  #   adjustfragments <- isolate(input$AdjustFragments)
  #   oldmod <- isolate(input$oldMod)
  #   newmod <- isolate(input$newMod)
  #   
  #   # Convert the .msp file to .db using the DBbuilder function
  #   db_path <- tempfile(fileext = ".db")
  #   DBbuilder(Library=Library(), FragmentationMode=FragmentationMode, MassAnalyzer=MassAnalyzer, CollisionEnergy=CollisionEnergy, CompoundClass=CompoundClassArg,
  #             Filter=Filter, DBoutput=db_path, topX=top, TMTPro = FALSE, cutoff=cutoff, massOffset=massOff, IonTypes=IonTypes, deltaFragment=adjustfragments, oldMod=oldmod, newMod=newmod)
  #   db_created(db_path)
  #   }
  #   # Read a random entry from the .db file
  #   conn <- dbConnect(RSQLite::SQLite(), db_created())
  #   random_entry <- dbGetQuery(conn, "SELECT * FROM CompoundTable INNER JOIN SpectrumTable ON CompoundTable.CompoundId = SpectrumTable.CompoundId ORDER BY RANDOM() LIMIT 1")
  #   dbDisconnect(conn)
  #   plotentry<<-data.frame(mass= 
  #                           lapply(random_entry$blobMass, function(x) {
  #     con <- rawConnection((x))
  #     mz <- readBin(con, double(), n = length(x) / 8)
  #     close(con)
  #     return(mz)
  #   })[[1]], intensity =
  #   
  #   intensity <<- lapply(random_entry$blobIntensity, function(x) {
  #     con <- rawConnection((x))
  #     intens <- readBin(con, double(), n = length(x) / 8)
  #     close(con)
  #     return(intens)
  #   })[[1]])
  #   
  #   return(plotentry)
  # 
  # }
  # 
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
         adjustfragments = input$AdjustFragments
         oldmod <- input$oldMod
         newmod <- input$newMod
         })
        DBbuilder(Library=Library(), FragmentationMode=FragmentationMode, MassAnalyzer=MassAnalyzer(), CollisionEnergy=CollisionEnergy, CompoundClass=CompoundClassArg,
                           Filter=Filter, DBoutput=x, topX=top, TMTPro = FALSE, cutoff=cutoff, massOffset=massOff, IonTypes=IonTypes,deltaFragment=adjustfragments, oldMod=oldmod, newMod=newmod)
      
   
    } )
  # output$previewPlot <- renderPlot({
  #   # Check if the 'Preview Entry' button has been clicked
  #   if (is.null(input$previewButton)) {
  #     return(NULL)
  #   }
  #   
  #   # Get the random entry
  #   random_entry <- get_random_entry()
  #   
  #   # Check if a random entry is available
  #   if (is.null(random_entry)) {
  #     return(NULL)
  #   }
  #   
  #   # Get the mass and intensity arrays from the random entry
  #    mass <- plotentry$mass
  #    intensity <- plotentry$intensity
  #   
  #   # Create the plot
  #   plot(plotentry$mass, plotentry$intensity,xlim=c(0,1000),ylim=c(0,1), type = "h", xlab = "Mass", ylab = "Intensity")
                                                                                      # "ModString:", random_entry$meta$ModString, "\n",
                                                                                      # "Sequence:", random_entry$meta$Sequence, "\n",
                                                                                      # "RetentionTime:", random_entry$meta$RetentionTime, "\n",
                                                                                      # "PrecursorMass:", random_entry$meta$PrecursorMass, "\n",
                                                                                      # "CollisionEnergy:", random_entry$meta$CollisionEnergy, "\n",
                                                                                      # "FragmentationMode:", random_entry$meta$FragmentationMode, "\n",
                                                                                      # "MassAnalyzer:", random_entry$meta$MassAnalyzer))
  #})
  
  

    
}

# Run the application 
shinyApp(ui = ui, server = server)
