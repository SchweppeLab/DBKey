## function for converting MSP files to DB files for RTLS analyses
# author: Chris McGann


#Takes txt file input from Prosit/SpectraST and extracts infomation needed to build .db
LibraryParser <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutof,massOffset) {
  
  PrositLib<- Library
  
  PrositLib<-stri_replace_all_fixed(PrositLib,"FullName: ", "")
  PrositLib<-stri_replace_all_fixed(PrositLib,"AvePrecursorMz:", "")
  
  #getting library entry names, should be in "PEPTIDEK/2" format
  Names<- PrositLib[stri_detect_fixed(PrositLib, "Name: ")]
  Names <- str_remove_all(Names, "Name: ")
  
  #Names without charges
  PeptideSequence <- gsub('.{2}$', '', Names)
  
  #Retrieve mod string for further processing 
  ModString<- str_extract(PrositLib,"ModString=[^=]+")
  ModString<- ModString[!is.na(ModString)]
  ModString<- str_extract(ModString,"//[^=]+")
  ModString<- gsub('.{6}$', '', ModString)
  ModString<- gsub('^.{2}', '', ModString)
  ModString<- gsub(' ', '', ModString)
  
  #testing until I figure out how to read unimod xml
  unimodTable <- read.csv("~/Repos/MSPtoDB/testMods.csv")
  ModString<-stri_replace_all_fixed(ModString, unimodTable$mod, as.character(unimodTable$massshift))
  
  
  #Getting charge state + last aa for adding mods
  Charge<-as.numeric(str_sub(Names,-1 ))
  LastAA<-str_sub(Names,-3,-3 )
  isK <- if(Source == "Prosit"){
    LastAA=="K"
  } else {
    LastAA=="]"
  }
  
  #Getting list of precursor masses, listed differently between Prosit/SpectraST
  PrecursorMasses <- if(Source == "Prosit"){
    PrositLib[stri_detect_fixed(PrositLib, "MW:")]
  } else {
    PrositLib[stri_detect_fixed(PrositLib, "PrecursorMZ:")]
  }
  PrecursorMasses <- gsub("[[:alpha:]]|:", "", PrecursorMasses)
  
  
  #If adding TMTPro, add mass to precursor mass
  PrecursorMasses<- if(TMTPro == TRUE) {
    as.numeric(PrecursorMasses)+304.207146/Charge
  } else  {PrecursorMasses}
  
  PrecursorMasses[isK] <- if(TMTPro == TRUE) {
    PrecursorMasses[isK]+304.207146/Charge[isK]
  } else {  PrecursorMasses[isK]}
 
  
   #Getting indicies of where peak lists are for further processing
  peakindexes<-if(Source=="Prosit") {
    which(stri_detect_fixed(PrositLib,"peaks:"))+1
  } else {
    which(stri_detect_fixed(PrositLib,"Peaks:"))+1
  }
  
  #I cant remember what problem I was trying to solve with this so I dont want to delete it yet
  nameindexes<-c(which(stri_detect_fixed(PrositLib,"Name: "))[-1]-1,length(PrositLib))
  cut<-min(length(peakindexes),length(nameindexes),length(PrecursorMasses))
  nameindexes<-nameindexes[1:cut]
  peakindexes<-peakindexes[1:cut]
  PrecursorMasses<-PrecursorMasses[1:cut]
  Names<-Names[1:cut]
  LastAA<-LastAA[1:cut]
  isK<-isK[1:cut]
  
  #retrieve peak lists from the indicies
  PeakLists<- (mapply(function(x, y) {PrositLib[x:y]}, x = peakindexes, y = nameindexes))
  PeakLists<-mapply(function(x) {str_split(PeakLists[[x]], "\\t", simplify = T)}, x = seq(from=1,to=length(PeakLists), by=1))
  names(PeakLists)<-PeptideSequence
  
  bYFunction<-function(x){
    pep<-unlist(str_split(x,""))
    b<- seq(1,length(pep), by=1)
    return(data.frame(aa=pep, b=paste0("b",b),y=paste0("y",y)))
  }
  
  
  OrganizePeaks<- function(x) {
    masses<-as.numeric(x[,1])
    int<-as.numeric(x[,2])
    annotations<-(x[,3])
    masses<-masses[!is.na(masses)]
    int<-int[!is.na(int)]
    peakNum <- min(topX, length(masses))
    
    dt<-data.table(masses,int, annotations)
    dt<-dt[dt$int>0,]
    maxPeak<-max(dt$int)
    dt<-dt[(dt$int/maxPeak)>cutoff,]
    dt<-setorder(dt, -int)
    dt <-dt[1:peakNum,]
    dt<-setorder(dt, masses)
    return(dt)
  }
  
  PeaksDT<- map(PeakLists, OrganizePeaks)  
  PeakAnnotations <- map(PeaksDT, function(x){(x[,3])}  ) 
  names(PeakAnnotations) <- ""
  PeakAnnotations <- map(PeakAnnotations, function(x){gsub("0.0ppm", "", x[[1]])})
  PeakAnnotations <- map(PeakAnnotations, function(x){gsub( '[^\\^|[:alnum:]]', "", x)})
 

  #takes peak lists, sorts, filters, adds TMTPro masses and reports blobs of masses
  #I should unify these functions 

#makes list of blobs for masses and intensities 
  blobMass<- as_blob(flatten(map(PeaksDT, function(x) {(as_blob(packBits(numToBits(x$masses))))})))
  names(blobMass)<- NULL
  
  blobInt<- as_blob(flatten(map(PeaksDT, function(x) {(as_blob(packBits(numToBits(x$int))))})))
  names(blobInt) <- NULL

  
  PeakAnnotationsCollapse<- map(PeakAnnotations, function(x) { paste(x, collapse = ";")})
  
  Tags<-paste0("mods:",ModString," ", "ions:", PeakAnnotationsCollapse )
  
  
#Adds lists to master list as it iterates through  
 # blobMassMaster <<- append(blobMassMaster, (blobMass))
 # blobIntMaster <<- append(blobIntMaster, (blobInt))
 # PrecursorMassesMaster <<-append(PrecursorMassesMaster, (PrecursorMasses))
 # NamesMaster <<- append(NamesMaster, (Names))
 # TagMaster <<- append(TagMaster, Tags)

 parallelTable<- data.table(blobMass=blobMass, blobInt=blobInt, 
                            PrecursorMasses=PrecursorMasses,Names=Names,Tags=Tags)
 return(parallelTable)
 
}


#Builds DB
DBbuilder<- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, DBoutput, Source, topX, cutoff) {

#Reads text file in using vroom  
# Library <-vroom("PrositTesting.msp", col_names = "Lib", delim = "\n")
# Library <- Library$Lib

Library <-fread(Library, header = F, sep  = "\n")
Library <- Library$V1

#gets poisition of where entires begin for splitting
#I'm going to fix this
AllNames<-grep("Name: ", Library, fixed = T)
AllNames<-AllNames[seq(1, length(AllNames), min(length(AllNames)-1,5000))]

#takes full library and converts it into a list of 5000 entry chunks
LibraryList<- Map(function(i,j) Library[i:j], AllNames, cumsum(diff(c(AllNames, length(Library)+1))))

blobMassMaster <<- list()
blobIntMaster <<- list()
PrecursorMassesMaster <<-list()
NamesMaster <<- list()
TagMaster <<- list()

n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#iterates through list building other lists
for(i in 1:length(LibraryList)) {
  LibraryParser(LibraryList[[i]], FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutoff," ")

}
resultsTable<-foreach(
  i = 1:length(LibraryList),
  .combine = 'c',
  .packages = c(
                "stringr",
                "stringi",
                "tidyverse",
                "readr",
                "data.table",
                "purrr",
                "blob")
   # That is the main point. Source your Function File here.

  
  ) %dopar% {
    source("~/Repos/MSPtoDB/lib/LibraryParser.R")
    LibraryParser(LibraryList[[i]], FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutoff," ")
}


parallel::stopCluster(cl = my.cluster)

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

}

#
#unimod <- read_xml("http://www.unimod.org/xml/unimod.xml")



