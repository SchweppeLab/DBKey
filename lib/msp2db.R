## function for converting MSP files to DB files for RTLS analyses
# author: Chris McGann

MSPtoDB <- function(Libraries, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, DBoutput, Source, topX, cutoff) {
  
  
  #Concatenate each of the libraries together
  outList<-NULL
  for(i in 1:length(Libraries[,1]))
  {
	outList[[i]] <- read_lines(Libraries[[i, 'datapath']])
  }
  
  #Flatten down to a single list
  PrositLib<-unlist(outList, recursive = FALSE)
  
  PrositLib<-stri_replace_all_fixed(PrositLib,"FullName: ", "")
  PrositLib<-stri_replace_all_fixed(PrositLib,"AvePrecursorMz:", "")
  
  
  
  Names<- PrositLib[stri_detect_fixed(PrositLib, "Name: ")]
  
  Names <- str_remove_all(Names, "Name: ")
  Charge<-as.numeric(str_sub(Names,-1 ))
  
  LastAA<-str_sub(Names,-3,-3 )
  
  
  isK <- if(Source == "Prosit"){
    LastAA=="K"
  } else {
    LastAA=="]"
  }
  
  PrecursorMasses <- if(Source == "Prosit"){
    PrositLib[stri_detect_fixed(PrositLib, "MW:")]
  } else {
    PrositLib[stri_detect_fixed(PrositLib, "PrecursorMZ:")]
  }
  PrecursorMasses <- gsub("[[:alpha:]]|:", "", PrecursorMasses)
  
  
  
  PrecursorMasses<- if(TMTPro == TRUE) {
    as.numeric(PrecursorMasses)+304.207146/Charge
  } else  {PrecursorMasses}
  
  PrecursorMasses[isK] <- if(TMTPro == TRUE) {
    PrecursorMasses[isK]+304.207146/Charge[isK]
  } else {  PrecursorMasses}
  
  peakindexes<-if(Source=="Prosit") {
    which(stri_detect_fixed(PrositLib,"peaks:"))+1
  } else {
    which(stri_detect_fixed(PrositLib,"Peaks:"))+1
  }
  
  nameindexes<-c(which(stri_detect_fixed(PrositLib,"Name: "))[-1]-1,length(PrositLib))
  
  
  
  PeakLists<- (mapply(function(x, y) {PrositLib[x:y]}, x = peakindexes, y = nameindexes))
  PeakLists<-mapply(function(x) {str_split(PeakLists[[x]], "\\t", simplify = T)}, x = seq(from=1,to=length(PeakLists), by=1))
  names(PeakLists)<-LastAA
  
  
  
  
  
  getMassBlob <- function(x) {
    masses<-as.numeric(x[[1]][,1])
    int<-as.numeric(x[1][[1]][,2])
    masses<-masses[!is.na(masses)]
    int<-int[!is.na(int)]
    
    dt<-data.table(masses,int)
    maxPeak<-max(dt$int)
    dt<-dt[(dt$int/maxPeak)>cutoff,]
    dt<-setorder(dt, -int)
    dt <-dt[1:topX,]
    dt<-setorder(dt, masses)
    
    blob<-(as_blob(packBits(numToBits(dt$masses))))
    return(blob)
  }
  getIntBlob <- function(x) {
    masses<-as.numeric(x[[1]][,1])
    int<-as.numeric(x[1][[1]][,2])
    masses<-masses[!is.na(masses)]
    int<-int[!is.na(int)]
    
    dt<-data.table(masses,int)
    maxPeak<-max(dt$int)
    dt<-dt[(dt$int/maxPeak)>cutoff,]
    dt<-setorder(dt, -int)
    dt <-dt[1:topX,]
    dt<-setorder(dt, masses)
    
    blob<-(as_blob(packBits(numToBits(dt$int))))
    return(blob)
  }
  
  
  x<-PeakLists[1]
  
  getMassBlobTMTpro <- function(x) {
    TMTPro<-304.207146
    isK<- names(x) == ""
    masses<-as.numeric(x[[1]][,1])
    int<-as.numeric(x[1][[1]][,2])
    masses<-masses[!is.na(masses)]
    int<-int[!is.na(int)]
    
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
    dt <-dt[1:topX,]
    dt<-setorder(dt, masses)

    
    Massblob<-as_blob(packBits(numToBits(dt$masses)))
    
    return(Massblob)
  }
  
  
  getIntBlobTMTpro <- function(x) {
    TMTPro<-304.207146
    isK<- names(x) == ""
    masses<-as.numeric(x[[1]][,1])
    int<-as.numeric(x[1][[1]][,2])
    
    frag<-(x[[1]][,3])
    frag<-frag[frag != ""]
    fragcharge<-str_extract_all(frag,"\\^.",simplify = T)
    fragcharge[fragcharge==""] <-1
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
    dt <-dt[1:topX,]
    dt<-setorder(dt, masses)
    
    
    
    intblob<-as_blob(packBits(numToBits(dt$int)))
    
    return(intblob)
  }
  
  
  blobMass<-if(TMTPro == TRUE) {
    as_blob(lmap(PeakLists, getMassBlobTMTpro))
  } else  {
    as_blob(lmap(PeakLists, getMassBlob))
  }
  
  blobInt<- if(TMTPro == TRUE) {
    as_blob(lmap(PeakLists, getIntBlobTMTpro))
  } else  {
    as_blob(lmap(PeakLists, getIntBlob))
  }
  
  
  
  
  
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
  
  conn4 <- dbConnect(SQLite(),DBoutput)
  
  dbWriteTable(conn4,"CompoundTable", CompoundTable, overwrite = T, field.types = c(CompoundId="INTEGER PRIMARY KEY", Synonyms="BLOB_TEXT", Structure="BLOB_TEXT"))
  
  dbWriteTable(conn4, "SpectrumTable", SpectrumTable, overwrite = T, field.types = c( SpectrumId = "INTEGER PRIMARY KEY", CompoundId ="INTEGER REFERENCES [CompoundTable]", blobMass="BLOB", blobIntensity="BLOB", RetentionTime = "DOUBLE", ScanNumber ="INTEGER", PrecursorMass="DOUBLE", NeutralMass="DOUBLE", blobAccuracy = "BLOB",blobFlags = "BLOB",blobResolution = "BLOB",blobNoises = "BLOB",blobTopPeaks = "BLOB", Version = "INTEGER"))
  
  dbWriteTable(conn4, "HeaderTable", HeaderTable, overwrite = T, field.types = c(version = "INTEGER NOT NULL DEFAULT 0", CreationDate="TEXT",LastModifiedDate="TEXT", Description="TEXT", Company="TEXT", ReadOnly = "BOOL", UserAccess="TEXT", PartialEdits = "BOOL"))
  
  dbWriteTable(conn4, "MaintenanceTable", MaintenanceTable, overwrite = T, field.types = c(CreationDate="TEXT",NoofCompoundsModified = "INTEGER", Description="TEXT"))
  
  dbDisconnect(conn4)
}
