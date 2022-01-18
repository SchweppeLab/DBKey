## function for converting MSP files to DB files for RTLS analyses
# author: Chris McGann


#Takes txt file input from Prosit/SpectraST and extracts infomation needed to build .db
LibraryParser <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutoff) {
  
  PrositLib<- Library
  
  PrositLib<-stri_replace_all_fixed(PrositLib,"FullName: ", "")
  PrositLib<-stri_replace_all_fixed(PrositLib,"AvePrecursorMz:", "")
  
  #getting library entry names, should be in "PEPTIDEK/2" format
  Names<- PrositLib[stri_detect_fixed(PrositLib, "Name: ")]
  Names <- str_remove_all(Names, "Name: ")
  
  #Names without charges
  NamesTag <- gsub('.{2}$', '', Names)
  
  #Retrieve mod string for further processing 
  ModString<- str_extract(PrositLib,"ModString=[^=]+")
  ModString<- ModString[!is.na(ModString)]
  ModString<- str_extract(ModString,"//[^=]+")
  ModString<- gsub('.{6}$', '', ModString)
  ModString<- gsub('^.{2}', '', ModString)
  
  
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
  names(PeakLists)<-LastAA
  
  
  #takes peak lists, sorts, filters, and reports blobs of masses
  getMassBlob <- function(x) {
    masses<-as.numeric(x[[1]][,1])
    int<-as.numeric(x[1][[1]][,2])
    annotations<-(x[1][[1]][,3])
    masses<-masses[!is.na(masses)]
    int<-int[!is.na(int)]
    peakNum <- min(topX, length(masses))
    
    dt<-data.table(masses,int)
    maxPeak<-max(dt$int)
    dt<-dt[(dt$int/maxPeak)>cutoff,]
    dt<-setorder(dt, -int)
    dt <-dt[1:peakNum,]
    dt<-setorder(dt, masses)   
    
    blob<-(as_blob(packBits(numToBits(dt$masses))))
    
    return(blob)
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
  PeakAnnotations <- lapply(PeakAnnotations, function(x){gsub('.{2}$', "", x)})
  PeakAnnotations <- lapply(PeakAnnotations, function(x){gsub('^.{1}', "", x)})
  
  #takes peak lists, sorts, filters, adds TMTPro masses and reports blobs of masses
  #I should unify these functions 
  getMassBlobTMTpro <- function(x) {
    TMTPro<-304.207146
    isK<- names(x) == ""
    masses<-as.numeric(x[[1]][,1])
    int<-as.numeric(x[1][[1]][,2])
    masses<-masses[!is.na(masses)]
    int<-int[!is.na(int)]
    peakNum <- min(topX, length(masses))
    
    frag<-(x[[1]][,3])
    frag<-frag[frag != ""]
    fragcharge<-str_extract_all(frag,"\\^.",simplify = T)
    fragcharge[fragcharge==""] <- 1
    fragcharge[!fragcharge==""]<-as.numeric(str_remove(fragcharge,"\\^"))
    if(is_empty(fragcharge)) {fragcharge=1}
    isBion<-stri_detect_fixed(frag,"b")
    
    dt<-data.table(masses,as.numeric(fragcharge),isBion,int)
    dt$masses[dt$isBion]<-dt$masses[dt$isBion]+TMTPro/dt$V2[dt$isBion]
    dt$masses[isK & !dt$isBion] <- dt$masses[!dt$isBion]+TMTPro/dt$V2[!dt$isBion]
    
    dt<-data.table(masses,int)
    maxPeak<-max(dt$int)
    dt<-dt[(dt$int/maxPeak)>cutoff,]
    dt<-setorder(dt, -int)
    dt <-dt[1:peakNum,]
    dt<-setorder(dt, masses)
    
    
    Massblob<-as_blob(packBits(numToBits(dt$masses)))
    
    return(Massblob)
  }
  
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
 blobMassMaster <<- append(blobMassMaster, (blobMass))
 blobIntMaster <<- append(blobIntMaster, (blobInt))
 PrecursorMassesMaster <<-append(PrecursorMassesMaster, (PrecursorMasses))
 NamesMaster <<- append(NamesMaster, (Names))
 TagMaster <<- append(TagMaster, Tags)

}


#Builds DB
DBbuilder<- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutoff) {

#Reads text file in using vroom  
Library <-vroom(Library, col_names = "Lib", delim = "\n")
Library <- Library$Lib

#gets poisition of where entires begin for splitting
#I'm going to fix this
AllNames<-grep("Name: ", Library, fixed = T)
AllNames<-AllNames[seq(1, length(AllNames), min(length(AllNames)-1,5000))]

#takes full library and converts it into a list of 5000 entry chunks
LibraryList<- Map(function(i,j) Library[i:j], AllNames, cumsum(diff(c(AllNames, length(Library)+1))))

blobMassMaster <- list()
blobIntMaster <- list()
PrecursorMassesMaster <-list()
NamesMaster <- list()
TagMaster <- list()


#iterates through list building other lists
for(i in 1:length(LibraryList)) {
  LibraryParser(LibraryList[[i]], FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutoff)
  
}


#Builds tables
MasterCompoundTable <- data.table(
  CompoundId = seq(1, length(blobIntMaster), by=1), 
  Formula = "", 
  Name = unlist(NamesMaster),
  Synonyms = "",
  Tag = unlist(TagMaster),
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


MasterSpectrumTable <- data.table(
  SpectrumId = seq(1, length(blobMassMaster), by=1), 
  CompoundId = seq(1, length(blobMassMaster), by=1),
  mzCloudURL = "",
  ScanFilter = "",
  RetentionTime = 0.0,
  ScanNumber = 0,#
  PrecursorMass = unlist(PrecursorMassesMaster),
  NeutralMass= 0,
  CollisionEnergy= CollisionEnergy,
  Polarity= "+",
  FragmentationMode= FragmentationMode,
  IonizationMode= "ESI",
  MassAnalyzer= MassAnalyzer,
  InstrumentName= "",
  InstrumentOperator= "",
  RawFileURL= "",
  blobMass=  I(blobMassMaster),
  blobIntensity= I(blobIntMaster),
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


DBbuilder(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutoff)
  
Library = "PrositTesting.msp"
FragmentationMode = "CID"
MassAnalyzer = "IT"
CollisionEnergy = "35"
Source = "Prosit"
TMTPro = FALSE
DBoutput = "TagTesting50klib.db"
topX = 50
cutoff = 0


