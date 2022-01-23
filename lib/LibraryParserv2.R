

LibraryParser <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, Filter, TMTPro, Source, topX, cutoff,massOffset) {
  library(stringr)                                                                                                                                                                                          
  library(stringi)                                                                                                                                                                                          
  library(data.table,warn.conflicts = FALSE )                                                                                                                                                                                       
  suppressPackageStartupMessages(library(purrr))                                                                                                                                                                                            
  library(blob)
  library(readr)

  nameindexes<-c(which(stri_detect_fixed(Library,"Name: ")))
  headerLength<- which(stri_detect_fixed(Library[1:100],"peaks:"))[1]-nameindexes[1]
  peakindexes<-nameindexes+headerLength
  HeaderLists<- unlist(mapply(function(x, y) {Library[x:y]}, x = nameindexes, y = peakindexes, SIMPLIFY = T))
  PeakLists<- mapply(function(x, y) {Library[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(Library)))
  
  HeaderLists<-gsub("FullName: ", "", HeaderLists, fixed = T)
  HeaderLists<-gsub("AvePrecursorMz:", "", HeaderLists, fixed = T)
  
  Names<- HeaderLists[stri_detect_fixed(HeaderLists, "Name: ")]
  Names <- str_remove_all(Names, "Name: ")
  NumPeaks<- HeaderLists[stri_detect_fixed(HeaderLists, "peaks: ")]
  NumPeaks <- gsub("[^[:digit:]]", "",  NumPeaks, perl = T)
  
  
  PeptideSequence <- gsub('.{2}$', '', Names, perl = T)
  
  #Retrieve mod string for further processing 
  ModString<- stri_extract_first_regex(HeaderLists,"ModString=[^=]+")
  ModString<- ModString[!is.na(ModString)]
  ModString<- str_extract(ModString,"//[^=]+")
  ModString<- gsub('.{6}$', '', ModString, perl = T)
  ModString<- gsub('^.{2}', '', ModString, perl = T)
  ModString<- gsub(' ', '', ModString, perl = T)

  unimodTable <- read.csv("~/Repos/MSPtoDB/testMods.csv")
  ModString<-stri_replace_all_fixed(ModString, unimodTable$mod, as.character(unimodTable$massshift), vectorize_all = F)
  

  Charge<-as.numeric(str_sub(Names,-1 ))
  LastAA<-str_sub(Names,-3,-3 )
  isK <- if(Source == "Prosit"){
    LastAA=="K"
  } else {
    LastAA=="]"
  }
  

  PrecursorMasses <- if(Source == "Prosit"){
    HeaderLists[stri_detect_fixed(HeaderLists, "MW:")]
  } else {
    HeaderLists[stri_detect_fixed(HeaderLists, "PrecursorMZ:")]
  }
  PrecursorMasses <- gsub("[[:alpha:]]|:", "", PrecursorMasses)
  
    PrecursorMasses<- if(TMTPro == TRUE) {
    as.numeric(PrecursorMasses)+304.207146/Charge
  } else  {PrecursorMasses}
  
  PrecursorMasses[isK] <- if(TMTPro == TRUE) {
    PrecursorMasses[isK]+304.207146/Charge[isK]
  } else {  PrecursorMasses[isK]}
  
rm(HeaderLists)
  PeakLists<- mapply(function(x, y) {Library[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(Library)),SIMPLIFY = T)
  PeakLists<-mapply(function(x) {str_split(PeakLists[[x]], "\\t",simplify = T)}, x = seq(from=1,to=length(PeakLists), by=1))
  PeakLists <- do.call("rbind", PeakLists)

  dt<-data.table(PeakLists)
  colnames(dt) <- c("masses", "int", "annotations")
  dt$masses<-as.numeric(dt$masses)
  dt$int<-as.numeric(dt$int)
  
  dt<-dt[dt$int>0,]
  rm(PeakLists)
  
  dt$annotations<- gsub( "[^//]+$","",dt$annotations, perl = T )
  dt$annotations<- gsub('[^\\^|[:alnum:]]', "", dt$annotations, perl=TRUE)
  
  PeakAnnotations<-mapply(function(x,y) {as.list(dt)[[3]][x:y]}, 
    x= c(0, cumsum(NumPeaks[-length(Names)])) ,y= cumsum(NumPeaks), SIMPLIFY = F)
  
  PeakAnnotations<-  lapply(PeakAnnotations, function(x) { paste(x, collapse = ";")})
  
  if(Filter == F){
  blobMass<-(as_blob(flatten(mapply(function(x,y) {
    as_blob(packBits(numToBits(as.list(dt)[[1]][x:y])))}, 
    x= c(0, cumsum(NumPeaks[-length(Names)])),y= cumsum(NumPeaks), SIMPLIFY = F))))
  
  
  blobInt<- (as_blob(flatten(mapply(function(x,y) {
    as_blob(packBits(numToBits(as.list(dt)[[2]][x:y])))}, 
    x= c(0, cumsum(NumPeaks[-length(Names)])),y= cumsum(NumPeaks), SIMPLIFY = F))))
  } else{
    PeakDT<-mapply(function(x,y) {dt[x:y]}, 
                            x= c(0, cumsum(NumPeaks[-length(Names)])) ,y= cumsum(NumPeaks), SIMPLIFY = F)
    OrganizePeaks<- function(x) {
      dt<-data.table(masses=x[,1],int=x[,2],annotations=x[,3])
      dt$masses<-as.numeric(substring(dt$masses,2))
      dt$int<-as.numeric(dt$int)
      
      dt<-dt[which(dt$int >0 ),]
      maxPeak<-max(dt$int)
      # dt<-dt[(dt$int/maxPeak)>cutoff,]
      if(topX < length(dt$masses)) {
        peakNum <- min(topX, length(dt$masses))
        dt<-setorder(dt, -int)
        dt <-dt[1:peakNum,]
      }
      #dt<-setorder(dt, masses)
      return(dt)
    }
    PeaksDT<- lapply(PeakLists, function(x){OrganizePeaks(x)}) 
    
    
  }
  
  
  Tags<-paste0("mods:",ModString," ", "ions:", PeakAnnotations )
 
   parallelTable<- data.table(blobMass=blobMass, blobInt=blobInt, 
                             PrecursorMasses=PrecursorMasses,Names=Names,Tags=Tags)

  return(parallelTable)
  
}

