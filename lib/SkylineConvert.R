library(RSQLite)
library(dplyr)
library(data.table)
OrganizePeaks<- function(x,topX,cutoff) {
  dt<-data.table(x)
  dt<-dt[which(dt$int >0 ),]

  if(topX < length(dt$masses)) {
    peakNum <- min(topX, length(dt$masses))
    dt$rank<- frank(dt$int)
    dt<-dt[dt$rank<=topX,]
    dt[,rank:=NULL]
  }
  if(cutoff>0) {
    maxPeak<-max(dt$int)
    dt<-dt[(dt$int/maxPeak)*100>cutoff,]
  }
  
  dt<-setkey(dt, masses)
  return(dt)
}
SkylineConvert<- function(x,CollisionEnergy,FragmentationMode,MassAnalyzer,topX,cutoff, Filter,massOffset, DBoutput) {


blib<-dbConnect(SQLite(),x)
mods<-dbReadTable(blib, "Modifications")
modsOutput <- mods %>% group_by(RefSpectraID) %>% summarise(ModString=paste0(mass, "@", position,";",collapse = ""))
RefSpectraPeaks<-dbReadTable(blib, "RefSpectraPeaks")
allmz<-lapply(RefSpectraPeaks$peakMZ, function(x) {readBin(memDecompress(unlist(x),'g'), double(), 200)})
allint<- lapply(RefSpectraPeaks$peakIntensity, function(x) {readBin(as.raw(unlist(x)),numeric(), n=200, size = 4)})
mzinttable <- lapply(seq(1,length(allmz), by = 1), function(x) {bind_cols(masses=allmz[[x]], int=allint[[x]])})

if(Filter) {
mzinttable<-lapply(seq(1,length(allmz), by =1), function(x) {OrganizePeaks(mzinttable[[x]],topX,cutoff)})
}



SpectraTable<-dbReadTable(blib, "RefSpectra")

allmzpacked<-lapply(mzinttable, function(x){(packBits(numToBits(unlist(x$masses))))})
allintpacked<-lapply(mzinttable, function(x){(packBits(numToBits(unlist(x$int))))})




if(length(massOffset) != 0){
  seqCharge <- data.frame(Sequence=SpectraTable$peptideSeq, charge =SpectraTable$precursorCharge, mZ= SpectraTable$precursorMZ )

  joined <-dplyr::left_join(seqCharge, read.csv(massOffset$datapath),by = "Sequence")
  
  joined$massOffsetTag<- ""
  joined$massOffsetTag[!is.na(joined$massOffset)] <- paste0("massOffset:",joined$massOffset[!is.na(joined$massOffset)])
  Tags<-paste0(joined$massOffsetTag," mods:",modsOutput$ModString," ", "ions:")
}
else {
  Tags<-paste0("mods:", modsOutput$ModString, " ", "ions:")

}



MasterCompoundTable <- data.table(
  CompoundId = seq(1, length(SpectraTable$id), by=1), 
  Formula = unlist(modsOutput$ModString), 
  Name = unlist(paste0(SpectraTable$peptideSeq, "/", SpectraTable$precursorCharge)),
  Synonyms = "",
  Tag = unlist(Tags),
  Sequence = SpectraTable$peptideSeq,
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
  SpectrumId = seq(1, length(SpectraTable$id), by=1), 
  CompoundId =seq(1, length(SpectraTable$id), by=1),
  mzCloudURL = "",
  ScanFilter = "",
  RetentionTime = SpectraTable$retentionTime,
  ScanNumber = 0,#
  PrecursorMass = SpectraTable$precursorMZ,
  NeutralMass= 0,
  CollisionEnergy= CollisionEnergy,
  Polarity= "+",
  FragmentationMode= FragmentationMode,
  IonizationMode= "ESI",
  MassAnalyzer= MassAnalyzer,
  InstrumentName= "",
  InstrumentOperator= "",
  RawFileURL= "",
  blobMass= I(allmzpacked), 
  blobIntensity= I(allintpacked), 
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

conn4 <- dbConnect(SQLite(),DBoutput)
dbWriteTable(conn4,"CompoundTable", MasterCompoundTable, overwrite = TRUE, append = F, field.types = c(CompoundId="INTEGER PRIMARY KEY", Synonyms="BLOB_TEXT", Structure="BLOB_TEXT"))
dbWriteTable(conn4, "SpectrumTable", MasterSpectrumTable, overwrite = TRUE, append = F, field.types = c( SpectrumId = "INTEGER PRIMARY KEY", CompoundId ="INTEGER REFERENCES [CompoundTable]", blobMass="BLOB", blobIntensity="BLOB", RetentionTime = "DOUBLE", ScanNumber ="INTEGER", PrecursorMass="DOUBLE", NeutralMass="DOUBLE", blobAccuracy = "BLOB",blobFlags = "BLOB",blobResolution = "BLOB",blobNoises = "BLOB",blobTopPeaks = "BLOB", Version = "INTEGER"))
dbWriteTable(conn4, "HeaderTable", HeaderTable, overwrite = TRUE, append = F, field.types = c(version = "INTEGER NOT NULL DEFAULT 0", CreationDate="TEXT",LastModifiedDate="TEXT", Description="TEXT", Company="TEXT", ReadOnly = "BOOL", UserAccess="TEXT", PartialEdits = "BOOL"))
dbWriteTable(conn4, "MaintenanceTable", MaintenanceTable, overwrite = TRUE, append = F, field.types = c(CreationDate="TEXT",NoofCompoundsModified = "INTEGER", Description="TEXT"))
dbDisconnect(conn4)

dbDisconnect(conn4)
}






