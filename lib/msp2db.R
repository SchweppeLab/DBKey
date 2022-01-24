## function for converting MSP files to DB files for RTLS analyses
# author: Chris McGann

#Builds DB
DBbuilder<- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro,Filter, DBoutput, topX, cutoff) {

  fileType<- str_extract(Library, "[^.]+$" )
if(fileType == "msp") {
  Source <- "Prosit"
} else if (fileType == "sptxt") {
  Source <- "SpectraST"
} else if( fileType == "blib") {
  Source <- "Skyline"
}
#Lazy load of library
LibraryRead <-vroom(Library, col_names = "Lib", delim = "\n")
LibraryRead<-LibraryRead$Lib

#List of interval to search for entry beginnings 
NamesX<-seq(200000,length(LibraryRead),by=200000)
NamesY<-seq(205000,length(LibraryRead),by=200000)

#Get indices of entries between 500 lines at evenly spaced intervals 
NamesList<-mapply(function(x, y) {(grep("Name: ", LibraryRead[x:y],fixed = T)+(x-1))[1]}, x = NamesX, y = NamesY)

#for function convienence
NamesListX<-c(0,NamesList)
NamesListY<-c(NamesList-1,length(LibraryRead))

#get rid of library
rm(LibraryRead)

#parallel with 8 cores. Load readr function and LibraryParser then read directly from file
# EightSnowParam<-paste0("8 core bplapply: ", system.time(bpmapply(function(x,y) {
#   source("~/Repos/MSPtoDB/lib/LibraryParser.R")
#   library(readr)
#   
#   Lib<-read_lines("myPrositLib.msp", skip=x-1, n_max=(y-x))
#   LibraryParser(Lib,"CID", "OT", "35", FALSE, "Prosit", 100, 0," ")},
#   x=NamesListX, y=NamesListY,
#          BPPARAM=SnowParam(workers = 8)))[3])
# 
# EightSnowParam<-paste0("8 core bplapply: ", system.time(bpmapply(function(x,y) {
#   source("~/Repos/MSPtoDB/lib/LibraryParserv2.R")
#   library(readr)
# 
#   Lib<-read_lines("myPrositLib.msp", skip=x-1, n_max=(y-x))
#   LibraryParser(Lib,"CID", "OT", "35", FALSE, "Prosit", 100, 0," ")},
#   x=NamesListX[1:100], y=NamesListY[1:100],
#   BPPARAM=SnowParam(workers = 8)))[3])

resultsTable<-bpmapply(function(x,y,z) {
  source("~/Repos/MSPtoDB/lib/LibraryParserv2.R")
  library(readr)
  
  Lib<-read_lines(z, skip=x-2, n_max=(y-x-1))
  LibraryParser(Lib,FragmentationMode, MassAnalyzer, CollisionEnergy, Filter, FALSE, Source, topX, cutoff," ")},
  x=NamesListX, y=NamesListY, z=Library,
  BPPARAM=SnowParam(workers = 10),SIMPLIFY = FALSE) %>% bind_rows


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
  CompoundClass= "",
  SmilesDescription= "",
  InChiKey= "")


MasterSpectrumTable <- data.table(
  SpectrumId = seq(1, length(resultsTable$Names), by=1), 
  CompoundId = seq(1, length(resultsTable$Names), by=1),
  mzCloudURL = "",
  ScanFilter = "",
  RetentionTime = 0.0,
  ScanNumber = 0,#
  PrecursorMass = resultsTable$PrecursorMasses,
  NeutralMass= 0,
  CollisionEnergy= CollisionEnergy,
  Polarity= "+",
  FragmentationMode= FragmentationMode,
  IonizationMode= "ESI",
  MassAnalyzer= MassAnalyzer,
  InstrumentName= "",
  InstrumentOperator= "",
  RawFileURL= "",
  blobMass= I(resultsTable$blobMass), #I(blobMassMaster),
  blobIntensity= I(resultsTable$blobInt), #I(blobIntMaster),
  blobAccuracy="",
  blobResolution="",
  blobNoises="",
  blobFlags="",
  blobTopPeaks="",
  Version="",
  CreationDate="",
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
  CreationDate =NA,
  NoofCompoundsModified=NA,
  Description=NA
)

##Writing new DB
conn4 <- dbConnect(SQLite(),DBoutput)
dbWriteTable(conn4,"CompoundTable", MasterCompoundTable, overwrite = T, append = F, field.types = c(CompoundId="INTEGER PRIMARY KEY", Synonyms="BLOB_TEXT", Structure="BLOB_TEXT"))
dbWriteTable(conn4, "SpectrumTable", MasterSpectrumTable, overwrite = T, append = F, field.types = c( SpectrumId = "INTEGER PRIMARY KEY", CompoundId ="INTEGER REFERENCES [CompoundTable]", blobMass="BLOB", blobIntensity="BLOB", RetentionTime = "DOUBLE", ScanNumber ="INTEGER", PrecursorMass="DOUBLE", NeutralMass="DOUBLE", blobAccuracy = "BLOB",blobFlags = "BLOB",blobResolution = "BLOB",blobNoises = "BLOB",blobTopPeaks = "BLOB", Version = "INTEGER"))
dbWriteTable(conn4, "HeaderTable", HeaderTable, overwrite = T, append = F, field.types = c(version = "INTEGER NOT NULL DEFAULT 0", CreationDate="TEXT",LastModifiedDate="TEXT", Description="TEXT", Company="TEXT", ReadOnly = "BOOL", UserAccess="TEXT", PartialEdits = "BOOL"))
dbWriteTable(conn4, "MaintenanceTable", MaintenanceTable, overwrite = T, append = F, field.types = c(CreationDate="TEXT",NoofCompoundsModified = "INTEGER", Description="TEXT"))
dbDisconnect(conn4)

dbDisconnect(conn4)
#return(conn4)
}

#
#unimod <- read_xml("http://www.unimod.org/xml/unimod.xml")



