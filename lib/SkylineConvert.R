SkylineConvert<- function(x, CollisionEnergy, MassAnalyzer, MassAnalyzer, DBoutput) {
  blib<-dbConnect(SQLite(), x)
  RefSpectra<-dbReadTable(blib, "RefSpectra")
  Modifications<-dbReadTable(blib, "Modifications")
  RefSpectraPeaks<-dbReadTable(blib, "RefSpectraPeaks")
  RefSpectraPeaks$numPeaks<-RefSpectra$numPeaks
  RT<-dbReadTable(blib, "RetentionTimes")
  
  RefSpectraPeaks$peakMZ <- lapply(RefSpectraPeaks$peakMZ, function(x) {as.raw(unlist(x))})
  RefSpectraPeaks$peakIntensity <-lapply(RefSpectraPeaks$peakMZ, function(x) {as.raw(unlist(x))})
  x<-RefSpectraPeaks[1,]
  
  ##From specL package
  convert_blib2psmInternal <- function(x,y,z){
    peakmZ<-(x)
    peakIntensity<-(y)
    numPeaks<-z
    
    mZ <- try(readBin(memDecompress(unlist(peakmZ),'g'), double(), numPeaks), TRUE)
    
    if (!is.numeric(mZ)|is_empty(mZ)){
      mZ <- try(readBin(as.raw(unlist(peakmZ)), double(),numPeaks), FALSE)
    }
    
    intensity <- try(readBin(memDecompress(unlist(peakIntensity),'g'), 
                             numeric(), n=numPeaks, size = 4), TRUE)
    
    
    if (!is.numeric(intensity) || length(intensity) != length(mZ)){
      intensity <- try(readBin(as.raw(unlist(peakIntensity)), 
                               numeric(), n=numPeaks, size = 4), FALSE)
      
    }  
    return(data.table(mZ=(as_blob(packBits(numToBits(mZ)))), int=(as_blob(packBits(numToBits(intensity))))))
  }
  zz<- mapply(function(x,y,z){convert_blib2psmInternal(x,y,z)}, x=RefSpectraPeaks[,2],
              y=RefSpectraPeaks[,3], RefSpectraPeaks[,4], SIMPLIFY = F)
  zz <- do.call("rbind", zz)
  
  Tags<-data.table(rep("mods:", length(RefSpectraPeaks$numPeaks)))
  addMods<-function(x) {
    ID<-x[[2]]
    pos<-x[3]
    mass<-x[4]
    y<-Tags[ID]
    Tags[ID]<<-paste0(y, mass, "@", pos," ")
  }
  for(i in 1:length(Modifications$id)) {
    addMods(Modifications[i,])
  }
  
  MasterCompoundTable <- data.table(
    CompoundId = RefSpectra$id,#seq(1, length(RefSpectra$id), by=1), 
    Formula = "", 
    Name = paste0(RefSpectra$peptideSeq,"/", RefSpectra$precursorCharge),
    Synonyms = "",
    Tag =unlist(Tags),
    Sequence =RefSpectra$peptideSeq,
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
  
  MasterSpectrumTable <- data.table(
    SpectrumId = RefSpectra$id,
    CompoundId =RefSpectra$id,
    mzCloudURL = NA,
    ScanFilter = NA,
    RetentionTime = RT$retentionTime,
    ScanNumber = NA,
    PrecursorMass = RefSpectra$precursorMZ,
    NeutralMass=NA,
    CollisionEnergy= CollisionEnergy,
    Polarity= "+",
    FragmentationMode= FragmentationMode,
    IonizationMode= NA,
    MassAnalyzer= MassAnalyzer,
    InstrumentName= NA,
    InstrumentOperator= NA,
    RawFileURL= NA,
    blobMass= I(zz$mZ), 
    blobIntensity= I(zz$int), 
    blobAccuracy=NA,
    blobResolution=NA,
    blobNoises=NA,
    blobFlags=NA,
    blobTopPeaks=NA,
    Version=NA,
    CreationDate=NA,
    Curator=NA,
    CurationType=NA,
    PrecursorIonType=NA,
    Acession=NA)
  
  ##Writing new DB
  conn4 <- dbConnect(SQLite(),"skylineTest1.db")
  dbWriteTable(conn4,"CompoundTable", MasterCompoundTable, overwrite = T, append = F, field.types = c(CompoundId="INTEGER PRIMARY KEY", Synonyms="BLOB_TEXT", Structure="BLOB_TEXT"))
  dbWriteTable(conn4, "SpectrumTable", MasterSpectrumTable, overwrite = T, append = F, field.types = c( SpectrumId = "INTEGER PRIMARY KEY", CompoundId ="INTEGER REFERENCES [CompoundTable]", blobMass="BLOB", blobIntensity="BLOB", RetentionTime = "DOUBLE", ScanNumber ="INTEGER", PrecursorMass="DOUBLE", NeutralMass="DOUBLE", blobAccuracy = "BLOB",blobFlags = "BLOB",blobResolution = "BLOB",blobNoises = "BLOB",blobTopPeaks = "BLOB", Version = "INTEGER"))
  dbWriteTable(conn4, "HeaderTable", HeaderTable, overwrite = T, append = F, field.types = c(version = "INTEGER NOT NULL DEFAULT 0", CreationDate="TEXT",LastModifiedDate="TEXT", Description="TEXT", Company="TEXT", ReadOnly = "BOOL", UserAccess="TEXT", PartialEdits = "BOOL"))
  dbWriteTable(conn4, "MaintenanceTable", MaintenanceTable, overwrite = T, append = F, field.types = c(CreationDate="TEXT",NoofCompoundsModified = "INTEGER", Description="TEXT"))
  
  
  dbDisconnect(conn4)
  
  
  
}


