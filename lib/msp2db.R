## function for converting MSP files to DB files for RTLS analyses
# author: Chris McGann

#Builds DB
library(BiocParallel)
library(RSQLite)
library(dplyr)
library(compiler)
library(vroom)
library(data.table)

# main function of shiny app
DBbuilder<- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy,
                     Filter, DBoutput, topX, cutoff, massOffset, IonTypes) {
# Get file type from input file
  fileType<- (Library$name[1])
  if(grepl("msp", fileType)) {
    fileType <- "Prosit"
  } else if(grepl("sptxt", fileType)) {
    fileType <- "SpectraST"
  } else if(grepl("blib", fileType)) {
    fileType<-"Skyline"
  }
  LibraryPath = (Library$datapath)  
 
  
  # Reading file in 
  if(fileType != "Skyline") {
   for (i in 1:length(Library$datapath))  {
    LibraryPath<-Library$datapath
    LibraryRead = vroom(LibraryPath, col_names = "Lib", delim = "\n",skip_empty_rows = FALSE)
   LibraryRead<-LibraryRead$Lib
  chunksize<-min(length(LibraryRead),500000)
  
  NamesX<-seq(chunksize,length(LibraryRead)-500,by=chunksize)
  NamesY<-seq(chunksize+500,length(LibraryRead),by=chunksize)
  
  #Get indices of entries between 500 lines at evenly spaced intervals 
  NamesList<-mapply(function(x, y) {(grep("^Name: ", LibraryRead[x:y],fixed = FALSE, perl = TRUE)+(x-1))[1]}, x = NamesX, y = NamesY)
#  NamesNamesList<-mapply(function(x, y) {LibraryRead[x:y][grepl("^Name: ", LibraryRead[x:y],fixed = FALSE, perl = TRUE)][[1]]}, x = NamesX, y = NamesY)
  
  #for function convienence
  NamesListX<-c(0,NamesList)
  NamesListY<-c(NamesList-1,length(LibraryRead))
  
  #get rid of library
 
  
  # chunklist<- mapply(function(x,y) {
  #   LibraryRead[x:y]
  # }, x=NamesListX, y=NamesListY)
   rm(LibraryRead)
  
  if(fileType=="SpectraST"){
    resultsTable<-bpmapply( function(x,y,z) {
      source("~/Repos/MSPtoDB/lib/SpXLibraryParser.R")
      Lib<-""
      Lib<-fread(z,skip = x-1, nrows = (y-x), blank.lines.skip=FALSE, sep = "\n",  header = FALSE )
      SpXLibraryParser(Library=Lib$V1, FragmentationMode="CID", MassAnalyzer="IT", 
                       CollisionEnergy="35", 
                       Filter=FALSE, TMTPro=FALSE, Source="SpectraST", topX=150, 
                       cutoff=0, IonTypes=NA)
    },
    x=NamesListX, y=NamesListY, z=LibraryPath,
    BPPARAM=SnowParam(workers = parallel::detectCores()-8),SIMPLIFY = FALSE) %>% bind_rows
  }

  # rtFix <- as.numeric(resultsTable$RetentionTime)
  # rtFix<- (rtFix/max(rtFix)*120)
  # resultsTable$RetentionTime <- as.character(rtFix)
   if(fileType=="Prosit"){
     resultsTable<-bpmapply(function(x,y,z) {
       source("~/Repos/MSPtoDB/lib/LibraryParserv2.R")
       Lib<-fread(z, skip=x-1, nrows=(y-x),strip.white = TRUE,header = FALSE, sep= "\n")
       LibraryParser(Library=Lib$V1, FragmentationMode=FragmentationMode, MassAnalyzer=MassAnalyzer,
                     CollisionEnergy=CollisionEnergy,
                     Filter=Filter, TMTPro=FALSE, Source=fileType, topX=topX,
                     cutoff=cutoff,massOffset, IonTypes=IonTypes)
     },
     x=NamesListX, y=NamesListY, z=LibraryPath,
     BPPARAM=SnowParam(workers = parallel::detectCores()-8),SIMPLIFY = FALSE) %>% bind_rows

   }
}
}

  
  
 if(fileType!="Skyline"){
 #Builds tables
    MasterCompoundTable <- data.table(
      CompoundId = seq(1, length(resultsTable$Names), by=1), 
      Formula = "", 
      Name = unlist(resultsTable$Names),
      Synonyms = "",
      Tag = unlist(resultsTable$Tags),
      Sequence = gsub('.{2}$', '',resultsTable$Names),
      CASId = "",
      ChemSpiderId= "",
      HMDBId= "",
      KEGGId= "",
      PubChemId= "",
      Structure= "",
      mzCloudId= as.integer(NA),
      CompoundClass= resultsTable$CompoundClass,
      SmilesDescription= "",
      InChiKey= "")
    
    
    MasterSpectrumTable <- data.table(
      SpectrumId = seq(1, length(resultsTable$Names), by=1), 
      CompoundId = seq(1, length(resultsTable$Names), by=1),
      mzCloudURL = "",
      ScanFilter = "",
      RetentionTime = resultsTable$RetentionTime,
      ScanNumber = 0,#
      PrecursorMass = resultsTable$PrecursorMasses,
      NeutralMass= 0,
      CollisionEnergy= resultsTable$CollisionEnergy,
      Polarity= "+",
      FragmentationMode= "CID",
      IonizationMode= "ESI",
      MassAnalyzer= "IT",
      InstrumentName= "",
      InstrumentOperator= "",
      RawFileURL= "",
      blobMass= I(resultsTable$blobMass), 
      blobIntensity= I(resultsTable$blobInt), 
      blobAccuracy="",
      blobResolution="",
      blobNoises="",
      blobFlags="",
      blobTopPeaks="",
      Version="",
      CreationDate=NA,
      Curator="",
      CurationType="",
      PrecursorIonType="",
      Acession="")
    
    HeaderTable <- data.frame(
      version = 5, 
      CreationDate = NA, 
      LastModifiedDate = NA,
      Description = NA,
      Company = NA,
      ReadOnly = NA,
      UserAccess = NA,
      PartialEdits= NA)
    
    MaintenanceTable <- data.frame(
      CreationDate = gsub("-"," ",Sys.Date()),
      NoofCompoundsModified=NA,
      Description=NA
    )
    
    ##Writing new DB
    conn4 <- dbConnect(SQLite(),DBoutput)
    dbWriteTable(conn4,"CompoundTable", MasterCompoundTable, overwrite = TRUE, append = F, field.types = c(CompoundId="INTEGER PRIMARY KEY", Synonyms="BLOB_TEXT", Structure="BLOB_TEXT"))
    dbWriteTable(conn4, "SpectrumTable", MasterSpectrumTable, overwrite = TRUE, append = F, field.types = c( SpectrumId = "INTEGER PRIMARY KEY", CompoundId ="INTEGER REFERENCES [CompoundTable]", blobMass="BLOB", blobIntensity="BLOB", RetentionTime = "DOUBLE", ScanNumber ="INTEGER", PrecursorMass="DOUBLE", NeutralMass="DOUBLE", blobAccuracy = "BLOB",blobFlags = "BLOB",blobResolution = "BLOB",blobNoises = "BLOB",blobTopPeaks = "BLOB", Version = "INTEGER"))
    dbWriteTable(conn4, "HeaderTable", HeaderTable, overwrite = TRUE, append = F, field.types = c(version = "INTEGER NOT NULL DEFAULT 0", CreationDate="TEXT",LastModifiedDate="TEXT", Description="TEXT", Company="TEXT", ReadOnly = "BOOL", UserAccess="TEXT", PartialEdits = "BOOL"))
    dbWriteTable(conn4, "MaintenanceTable", MaintenanceTable, overwrite = TRUE, append = F, field.types = c(CreationDate="TEXT",NoofCompoundsModified = "INTEGER", Description="TEXT"))
    dbDisconnect(conn4)
    
    dbDisconnect(conn4)
    #return(conn4)
  }
   if( fileType == "Skyline") {
    source("~/Repos/MSPtoDB/lib/SkylineConvert.R")
    
    SkylineConvert(LibraryPath, CollisionEnergy=CollisionEnergy, FragmentationMode=FragmentationMode, 
                  MassAnalyzer=MassAnalyzer, DBoutput=DBoutput) 
    
  }
}

# outList<-NULL
#  for(i in 1:length(Library$name))
#  {
#    outList[i] <- vroom(Library$datapath[i], col_names = "Lib", delim = "\n")
#  }
#  LibraryRead<-unlist(outList, recursive = FALSE)
# 
#     
# if(length(Library$name==1)){
#   LibraryPath = Library$datapath
# } else {
#    dir.create("Temp")
#    fwrite(LibraryRead, "CombinedLibTemp.msp")
#    LibraryPath = "Temp/CombinedLibTemp.msp"
#  }
#   
#  
#   fileType<- substring(Library$name, nchar(Library$name)-3)
# if(fileType == "msp") {
#   Source <- "Prosit"
# } else if (fileType == "ptxt") {
#   Source <- "SpectraST"
# } else if( fileType == "blib") {
#   Source <- "Skyline"
# }

