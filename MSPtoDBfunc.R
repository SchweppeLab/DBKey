#Packages
```{r}
library(RSQLite)
library(stringr)
library(Rcpp)
library(stringi)
library(dplyr)
library(readr)
library(data.table)
library(purrr)
library(blob)
library(OrgMassSpecR)

```


#auxiliary functions 
```{r}
#Takes list of Masses/Int/annotations returns blob of masses
getMassBlob <- function(x) {
  y<-(x[[1]][,1])
  blob<-as_blob(packBits(numToBits(y)))
  return(blob)
}

#Takes list of Masses/Int/annotations returns blob of masses with TMTPro added
#Doesnt work with missed cleavage. 
getMassBlobTMTpro <- function(x) {
  TMTPro<-304.207146
  isK<- names(x) == "K"
  masses<-as.numeric(x[[1]][,1])
  frag<-(x[[1]][,3])
  fragcharge<-str_extract_all(frag,"\\^.",simplify = T)
  fragcharge[fragcharge==""] <-1
  fragcharge[!fragcharge==""]<-as.numeric(str_remove(fragcharge,"\\^"))
  isBion<-stri_detect_fixed(frag,"b")
  
  dt<-data.table(masses,as.numeric(fragcharge),isBion)
  dt$masses[dt$isBion]<-dt$masses[dt$isBion]+TMTPro/dt$V2[dt$isBion]
  dt$masses[isK & !dt$isBion] <- dt$masses[!dt$isBion]+TMTPro/dt$V2[!dt$isBion]
  
  
  blob<-as_blob(packBits(numToBits(dt$masses)))
  return(blob)
}
#Takes list of Masses/Int/annotations returns blob of intensities
getIntBlob <- function(x) {
  y<-(x[[1]][,2])
  blob<-as_blob(packBits(numToBits(y)))
  return(blob)
}

```


#Read .msp and build tables 
```{r}
#Calls function. Library->.msp file path, CID/HCD, IT/FT, integer, boolean, name of .db output
#example: MSPtoDB("test.msp", "CID", "IT", 35, TRUE, "testingfun.db", SpectraST/Prosit )


MSPtoDB <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, DBoutput, Source) {
  #reads msp
  PrositLib<- read_lines(Library)
  
  # Pulls names, charge and last AA from file
  Names<- PrositLib[stri_detect_fixed(PrositLib, "Name: ")]
  Names <- str_remove_all(Names, "Name: ")
  Charge<-as.numeric(str_sub(Names,-1 ))
  LastAA<-str_sub(Names,-3,-3 )
  
  #Pulls list of precursor masses from file
  PrecursorMasses <- if(Source == "Prosit"){
    PrositLib[stri_detect_fixed(PrositLib, "MW:")]
  }
  else {
    PrositLib[stri_detect_fixed(PrositLib, "PrecursorMZ:")]
  }
  PrecursorMasses <- gsub("[[:alpha:]]|:", "", PrecursorMasses)
  
  #finds the index of the beginning and end of all peak lists
  peakindexes<-if(Source == "Prosit") {
    which(stri_detect_fixed(PrositLib,"peaks:"))+1}
  else {
    which(stri_detect_fixed(PrositLib,"Peaks:"))+1
  }
  nameindexes<-c(which(stri_detect_fixed(PrositLib,"Name: "))[-1]-1,length(PrositLib))
  
  #Subsets the file to make a list of all mass/intensity/annotation strings
  PeakLists<- (mapply(function(x, y) {PrositLib[x:y]}, x = peakindexes, y = nameindexes))
  
  #Splits the strings into matrix, adds last amino acid for later use
  PeakLists<-mapply(function(x) {str_split(PeakLists[[x]], "\\t", simplify = T)}, x = seq(from=1,to=length(PeakLists), by=1))
  names(PeakLists)<-LastAA
  
  #Removes things that are no longer needed
  PrositLib<-NULL
  peakindexes<-NULL
  nameindexes<-NULL
  
  
  #makes blobs of peaks and ints, adds TMT to frag if applicable
  blobMass<-if(TMTPro=="TRUE"){ 
    (as_blob(lmap(PeakLists, getMassBlobTMTpro)))}
  else {
    (as_blob(lmap(PeakLists, getMassBlob)))
  }
  
  blobInt<- as_blob(lmap(PeakLists, getIntBlob))
  
  #removes peaklists
  PeakLists<-NULL
  
  #Create compound table
  CompoundTable <- data.table(
    CompoundId = as.integer(seq(1, length(Names), by=1)), 
    Formula = "", 
    Name = Names,
    Synonyms = "",
    Tag = "",
    Sequence = "",
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
  
  #create spectrum table
  SpectrumTable <- data.table(
    SpectrumId = seq(1, length(Names), by=1), 
    CompoundId = seq(1, length(Names), by=1), 
    mzCloudURL = "",
    ScanFilter = "",
    RetentionTime = 0.0,
    ScanNumber = 0,#
    PrecursorMass = PrecursorMasses,
    NeutralMass= 0,
    CollisionEnergy= CollisionEnergy,
    Polarity= "+",
    FragmentationMode= FragmentationMode,
    IonizationMode= "ESI",
    MassAnalyzer= MassAnalyzer,
    InstrumentName= "",
    InstrumentOperator= "",
    RawFileURL= "",
    blobMass=  I(blobMass),
    blobIntensity= I(blobInt),
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
  
  #creates HeaderTable
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
  
  #creates .db
  filename<-as.character(DBoutput)
  
  conn3 <- dbConnect(SQLite(), dbname = DBoutput)
  
  #writes all tables to db. Overwrites if same db used.
  dbWriteTable(conn3,"CompoundTable", CompoundTable, overwrite = T, field.types = c(CompoundId="INTEGER PRIMARY KEY", Synonyms="BLOB_TEXT", Structure="BLOB_TEXT"))
  dbWriteTable(conn3, "SpectrumTable", SpectrumTable, overwrite = T, field.types = c( SpectrumId = "INTEGER PRIMARY KEY", CompoundId ="INTEGER REFERENCES [CompoundTable]", blobMass="BLOB", blobIntensity="BLOB", RetentionTime = "DOUBLE", ScanNumber ="INTEGER", PrecursorMass="DOUBLE", NeutralMass="DOUBLE", blobAccuracy = "BLOB",blobFlags = "BLOB",blobResolution = "BLOB",blobNoises = "BLOB",blobTopPeaks = "BLOB", Version = "INTEGER"))
  dbWriteTable(conn3, "HeaderTable", HeaderTable, overwrite = T, field.types = c(version = "INTEGER NOT NULL DEFAULT 0", CreationDate="TEXT",LastModifiedDate="TEXT", Description="TEXT", Company="TEXT", ReadOnly = "BOOL", UserAccess="TEXT", PartialEdits = "BOOL"))
  dbWriteTable(conn3, "MaintenanceTable", MaintenanceTable, overwrite = T, field.types = c(CreationDate="TEXT",NoofCompoundsModified = "INTEGER", Description="TEXT"))
  return(conn3)
  #disconnects database
  dbDisconnect(conn3)
}


