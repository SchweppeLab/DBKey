SkylineConvert<- function(x) {
  blib<-dbConnect(SQLite(), x)
  RefSpectra<-dbReadTable(blib, "RefSpectra")
  Modifications<-dbReadTable(blib, "Modifications")
  RefSpectraPeaks<-dbReadTable(blib, "RefSpectraPeaks")
  RetentionTimes<-dbReadTable(blib, "RetentionTimes")
  
  hm<-readBin(RefSpectraPeaks$peakIntensity[[4]], what='double', n = 100, size=4)
  hmm<-as_blob(as_blob(packBits(numToBits(hm))))

  MasterCompoundTable <- data.table(
    CompoundId = c(1), 
    Formula = "", 
    Name = RefSpectra$peptideSeq[4],
    Synonyms = "",
    Tag = "",
    Sequence = RefSpectra$peptideSeq[4],
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
    SpectrumId = c(1), 
    CompoundId =c(1),
    mzCloudURL = "",
    ScanFilter = "",
    RetentionTime = 0.0,
    ScanNumber = 0,#
    PrecursorMass = RefSpectra$precursorMZ[4],
    NeutralMass= 0,
    CollisionEnergy= CollisionEnergy,
    Polarity= "+",
    FragmentationMode= FragmentationMode,
    IonizationMode= "ESI",
    MassAnalyzer= MassAnalyzer,
    InstrumentName= "",
    InstrumentOperator= "",
    RawFileURL= "",
    blobMass=I(RefSpectraPeaks$peakMZ[4]), #I(blobMassMaster),
    blobIntensity= I(hmm), #I(blobIntMaster),
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
  conn4 <- dbConnect(SQLite(),"skylineTest.db")
  dbWriteTable(conn4,"CompoundTable", MasterCompoundTable, overwrite = T, append = F, field.types = c(CompoundId="INTEGER PRIMARY KEY", Synonyms="BLOB_TEXT", Structure="BLOB_TEXT"))
  dbWriteTable(conn4, "SpectrumTable", MasterSpectrumTable, overwrite = T, append = F, field.types = c( SpectrumId = "INTEGER PRIMARY KEY", CompoundId ="INTEGER REFERENCES [CompoundTable]", blobMass="BLOB", blobIntensity="BLOB", RetentionTime = "DOUBLE", ScanNumber ="INTEGER", PrecursorMass="DOUBLE", NeutralMass="DOUBLE", blobAccuracy = "BLOB",blobFlags = "BLOB",blobResolution = "BLOB",blobNoises = "BLOB",blobTopPeaks = "BLOB", Version = "INTEGER"))
  dbWriteTable(conn4, "HeaderTable", HeaderTable, overwrite = T, append = F, field.types = c(version = "INTEGER NOT NULL DEFAULT 0", CreationDate="TEXT",LastModifiedDate="TEXT", Description="TEXT", Company="TEXT", ReadOnly = "BOOL", UserAccess="TEXT", PartialEdits = "BOOL"))
  dbWriteTable(conn4, "MaintenanceTable", MaintenanceTable, overwrite = T, append = F, field.types = c(CreationDate="TEXT",NoofCompoundsModified = "INTEGER", Description="TEXT"))
  dbDisconnect(conn4)
  
  dbDisconnect(conn4)
  
  
  
  
  
  
}
