
library(stringr)                                                                                                                                                                                          
library(stringi)                                                                                                                                                                                          
library(data.table,warn.conflicts = FALSE )                                                                                                                                                                                       
suppressPackageStartupMessages(library(purrr))                                                                                                                                                                                            
library(blob)
library(readr)
library(compiler)
OrganizePeaks<- function(x,topX,cutoff,IonTypes) {
  
  dt<-x[which(x$int >0 ),]
  
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
  if(!is.null(IonTypes)) {
    dt<-dt[stri_detect_fixed(x$annotations,IonTypes),]
  }
  dt<-setkey(dt, masses)
  return(dt)
}

OrderPeaks<- function(x) {
  
  dt<-setkey(x, masses)
  return(dt)
}

blobMassFunction<- function(x){
  as_blob(packBits(numToBits(unlist(x[,1]))))}
blobIntFunction<- function(x){
  as_blob(packBits(numToBits(unlist(x[,2]))))}


OrgPeakCmp<-cmpfun(OrganizePeaks)
OrdPeakCmp<-cmpfun(OrderPeaks)
blobMassFunctionCmp<- cmpfun(blobMassFunction)
blobIntFunctionCmp<- cmpfun(blobIntFunction)






LibraryParser <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, 
                          Filter=F, TMTPro=F, Source, topX=0, cutoff=0,massOffset, IonTypes=NA) {

  
  nameindexes<-c(which(stri_detect_fixed(Library,"Name: ")))
  headerLength<- which(stri_detect_regex(Library[1:100],"(?i)peaks:"))[1]-nameindexes[1]
  peakindexes<-nameindexes+headerLength
  HeaderLists<- unlist(mapply(function(x, y) {Library[x:y]}, x = nameindexes, y = peakindexes, SIMPLIFY = T))
  PeakLists<- mapply(function(x, y) {Library[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(Library)))
  
  HeaderLists<-gsub("FullName: ", "", HeaderLists, fixed = T)
  HeaderLists<-gsub("AvePrecursorMz:", "", HeaderLists, fixed = T)
  
    
    getFrag<- function(x){
     match<- unique(stri_extract_all_fixed(x, c("CID","HCD"), simplify = T, omit_no_match = T))
     if(length(match)==2){
       return(as.character(match[1]))
     } else {
       stop("Unable to determine Fragmentation method from file, please specify")
     }
    }
    
    if(FragmentationMode== "Read From file"){
      FragmentationMode=getFrag(HeaderLists[1:100])
    }else{
      FragmentationMode=FragmentationMode }
    
    
    
    getCE<- function(x){
      
      CE<-x[stri_detect_fixed(x, "Collision")]
      CE<-gsub("_", "", CE)
      CE<- stri_extract_first_regex(CE,"Collisionenergy=[^/d]{2,4}")
      CE<- gsub("Collisionenergy=", "", CE)
      if(length(CE)==0) {
        stop("Unable to determine Collision Energies from file, please specify")
      } else{
        return(CE)
      }
    }
      
    
      CollisionEnergy<-  if(CollisionEnergy== "Read from file"){
      getCE(HeaderLists)
    } else{
      CollisionEnergy }
  
  Names<- HeaderLists[stri_detect_fixed(HeaderLists, "Name: ")]
  Names <- str_remove_all(Names, "Name: ")
  NumPeaks<-(c(nameindexes[-1], length(Library)+1) - peakindexes)-1

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
  # LastAA<-str_sub(Names,-3,-3 )
  # isK <- if(Source == "Prosit"){
  #   LastAA=="K"
  # } else {
  #   LastAA=="]"
  # }
  # 

  PrecursorMasses <- if(Source == "Prosit"){
    HeaderLists[stri_detect_fixed(HeaderLists, "MW:")]
  } else {
    HeaderLists[stri_detect_fixed(HeaderLists, "PrecursorMZ:")]
  }
  PrecursorMasses <- gsub("[[:alpha:]]|:", "", PrecursorMasses)
  
  #   PrecursorMasses<- if(TMTPro == TRUE) {
  #   as.numeric(PrecursorMasses)+304.207146/Charge
  # } else  {PrecursorMasses}
  # 
  # PrecursorMasses[isK] <- if(TMTPro == TRUE) {
  #   PrecursorMasses[isK]+304.207146/Charge[isK]
  # } else {  PrecursorMasses[isK]}
  
rm(HeaderLists)
  
  PeakLists<- mapply(function(x, y) {Library[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(Library)),SIMPLIFY = T)
  PeakLists<-mapply(function(x) {str_split(PeakLists[[x]], "\\t",simplify = T)}, x = seq(from=1,to=length(PeakLists), by=1))
  PeakLists <- do.call("rbind", PeakLists)

  dt<-data.table(PeakLists)
  colnames(dt) <- c("masses", "int", "annotations")
  dt$masses<-as.numeric(dt$masses)
  dt$int<-as.numeric(dt$int)
  
   UnSorted<- any(unlist(mapply(function(x,y) {is.unsorted(dt[[1]][x:y])}, 
          x= c(0, cumsum(NumPeaks[1:4]))+1 ,y= cumsum(NumPeaks[1:5]), SIMPLIFY = F)))
  
  rm(PeakLists)
  
  dt$annotations<- gsub( "[^//]+$","",dt$annotations, perl = T )
  dt$annotations<- gsub('[^\\^|[:alnum:]]', "", dt$annotations, perl=TRUE)
  

  if(Filter == F & !UnSorted){

    
   blobMass<-(as_blob(flatten(mapply(function(x,y,z){as_blob(packBits(numToBits(as.list(dt)[[1]][x:y])))},
     x=c(0, cumsum(NumPeaks[-length(Names)]))+1,y=cumsum(NumPeaks),SIMPLIFY = F))))
  
    blobInt<-(as_blob(flatten(mapply(function(x,y,z){as_blob(packBits(numToBits(as.list(dt)[[2]][x:y])))},
               x=c(0, cumsum(NumPeaks[-length(Names)]))+1,y=cumsum(NumPeaks),SIMPLIFY = F))))
    
  }     else if (Filter == T) {
    PeakDT<-mapply(function(x,y) {dt[x:y]},
          x= c(0, cumsum(NumPeaks[-length(Names)]))+1 ,y= cumsum(NumPeaks), SIMPLIFY = F)

    PeakDTOrganize<- lapply(PeakDT, function(x){OrgPeakCmp(x,topX,cutoff,IonTypes)})
    rm(PeakDT)

    blobMass<-(as_blob(flatten(mapply(blobMassFunctionCmp,PeakDTOrganize, SIMPLIFY = F))))

    blobInt<-(as_blob(flatten(mapply(blobIntFunctionCmp,PeakDTOrganize, SIMPLIFY = F))))
    
  } else {
    PeakDT<-mapply(function(x,y) {dt[x:y]},
                   x= c(0, cumsum(NumPeaks[-length(Names)]))+1 ,y= cumsum(NumPeaks), SIMPLIFY = F)

    PeakDTOrganize<- lapply(PeakDT, function(x){OrdPeakCmp(x)})
    rm(PeakDT)

    blobMass<-(as_blob(flatten(mapply(blobMassFunctionCmp,PeakDTOrganize, SIMPLIFY = F))))

    blobInt<-(as_blob(flatten(mapply(blobIntFunctionCmp,PeakDTOrganize, SIMPLIFY = F))))

  }
  
  PeakAnnotations<-mapply(function(x,y) {as.list(dt)[[3]][x:y]},
                          x= c(0, cumsum(NumPeaks[-length(Names)])+1) ,y= cumsum(NumPeaks), SIMPLIFY = F)

  PeakAnnotations<-  lapply(PeakAnnotations, function(x) { paste(x, collapse = ";")})
  
  
  if(massOffset==0){
    Tags<-paste0("mods:",ModString," ", "ions:", PeakAnnotations)
  } else {
    Tags<-paste0("MassOffset: ",massOffset," mods:",ModString," ", "ions:", PeakAnnotations )
  }

     parallelTable<- data.table(blobMass=blobMass, blobInt=blobInt, 
                             PrecursorMasses=PrecursorMasses,Names=Names,Tags=Tags, FragmentationMode=FragmentationMode,
                             CollisionEnergy=CollisionEnergy)

  return(parallelTable)
  
}

LibraryParserCmp<-cmpfun(LibraryParser)

SpXLibraryParser <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, 
                             Filter=F, TMTPro=F, Source, topX=0, cutoff=0,massOffset, IonTypes=NA) {
  
  firstName<-grep("Name", Library[1:100])[1]
  Library<-Library[firstName:length(Library)]
  
  
  nameindexes<-c(1,which(stri_detect_fixed(Library,"\nName: ")))
  headerLength<- which(stri_detect_regex(Library[1:100],"(?i)peaks:"))[1]-nameindexes[1]
  peakindexes<-nameindexes+headerLength
  HeaderLists<- unlist(mapply(function(x, y) {Library[x:y]}, x = nameindexes, y = peakindexes, SIMPLIFY = T))
  PeakLists<- mapply(function(x, y) {Library[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(Library)))
  
  
  
  getFrag<- function(x){
    match<- unique(stri_extract_all_fixed(x, c("CID","HCD"), simplify = T, omit_no_match = T))
    if(length(match)==2){
      return(as.character(match[1]))
    } else {
      stop("Unable to determine Fragmentation method from file, please specify")
    }
  }
  
  if(FragmentationMode== "Read From file"){
    getFrag(HeaderLists[1:100])
  }else{
    FragmentationMods=FragmentationMode }
  
  
  
  getCE<- function(x){
    
    CE<-x[stri_detect_fixed(x, "Collision")]
    CE<-gsub("_", "", CE)
    CE<- stri_extract_first_regex(CE,"Collisionenergy=[^/d]{2,4}")
    CE<- gsub("Collisionenergy=", "", CE)
    if(length(CE)==0) {
      stop("Unable to determine Collision Energies from file, please specify")
    } else{
      return(CE)
    }
  }
  
  
  CollisionEnergy<-  if(CollisionEnergy== "Read from file"){
    getCE(HeaderLists)
  } else{
    CollisionEnergy }
  
  HeaderLists<-gsub("FullName: ", "", HeaderLists, fixed = T)
  HeaderLists<-gsub("AvePrecursorMz:", "", HeaderLists, fixed = T)
  
  
  Names<- HeaderLists[stri_detect_fixed(HeaderLists, "Name: ")]
  Names <- str_remove_all(Names, "\nName: ")
  Names[1] <- str_remove_all(Names[1], "Name: ")
  NumPeaks<-(c(nameindexes[-1], length(Library)+1) - peakindexes)-1

  
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
  # LastAA<-str_sub(Names,-3,-3 )
  # isK <- if(Source == "Prosit"){
  #   LastAA=="K"
  # } else {
  #   LastAA=="]"
  # }
  # 
  
  PrecursorMasses <- if(Source == "Prosit"){
    HeaderLists[stri_detect_fixed(HeaderLists, "MW:")]
  } else {
    HeaderLists[stri_detect_fixed(HeaderLists, "PrecursorMZ:")]
  }
  PrecursorMasses <- gsub("[[:alpha:]]|:", "", PrecursorMasses)
  
  #   PrecursorMasses<- if(TMTPro == TRUE) {
  #   as.numeric(PrecursorMasses)+304.207146/Charge
  # } else  {PrecursorMasses}
  # 
  # PrecursorMasses[isK] <- if(TMTPro == TRUE) {
  #   PrecursorMasses[isK]+304.207146/Charge[isK]
  # } else {  PrecursorMasses[isK]}
  
  rm(HeaderLists)
  
  PeakLists<- mapply(function(x, y) {Library[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(Library)),SIMPLIFY = T)
  PeakLists<-mapply(function(x) {str_split(PeakLists[[x]], "\\t",simplify = T)}, x = seq(from=1,to=length(PeakLists), by=1))
  PeakLists <- do.call("rbind", PeakLists)
  
  dt<-data.table(PeakLists)
  dt<-dt[,1:3]
  colnames(dt) <- c("masses", "int", "annotations")
  dt$masses<-as.numeric(dt$masses)
  dt$int<-as.numeric(dt$int)
  
  UnSorted<- any(unlist(mapply(function(x,y) {is.unsorted(dt[[1]][x:y])}, 
                               x= c(0, cumsum(NumPeaks[1:4]))+1 ,y= cumsum(NumPeaks[1:5]), SIMPLIFY = F)))
  
  rm(PeakLists)
  
  dt$annotations<- gsub( "[^//]+$","",dt$annotations, perl = T )
  dt$annotations<- gsub('[^\\^|[:alnum:]]', "", dt$annotations, perl=TRUE)
  
  
  # UnSorted<-F
  if(Filter == F & !UnSorted){
    
    
    blobMass<-(as_blob(flatten(mapply(function(x,y,z){as_blob(packBits(numToBits(as.list(dt)[[1]][x:y])))},
                                      x=c(0, cumsum(NumPeaks[-length(Names)]))+1,y=cumsum(NumPeaks),SIMPLIFY = F))))
    
    blobInt<-(as_blob(flatten(mapply(function(x,y,z){as_blob(packBits(numToBits(as.list(dt)[[2]][x:y])))},
                                     x=c(0, cumsum(NumPeaks[-length(Names)]))+1,y=cumsum(NumPeaks),SIMPLIFY = F))))
    
  }     else if (Filter == T) {
    PeakDT<-mapply(function(x,y) {dt[x:y]},
                   x= c(0, cumsum(NumPeaks[-length(Names)]))+1 ,y= cumsum(NumPeaks), SIMPLIFY = F)
    
    PeakDTOrganize<- lapply(PeakDT, function(x){OrgPeakCmp(x,topX,cutoff,IonTypes)})
    rm(PeakDT)
    
    blobMass<-(as_blob(flatten(mapply(blobMassFunctionCmp,PeakDTOrganize, SIMPLIFY = F))))
    
    blobInt<-(as_blob(flatten(mapply(blobIntFunctionCmp,PeakDTOrganize, SIMPLIFY = F))))
    
  } else {
    PeakDT<-mapply(function(x,y) {dt[x:y]},
                   x= c(0, cumsum(NumPeaks[-length(Names)])) ,y= cumsum(NumPeaks), SIMPLIFY = F)
    
    PeakDTOrganize<- lapply(PeakDT, function(x){OrdPeakCmp(x)})
    rm(PeakDT)
    
    blobMass<-(as_blob(flatten(mapply(blobMassFunctionCmp,PeakDTOrganize, SIMPLIFY = F))))
    
    blobInt<-(as_blob(flatten(mapply(blobIntFunctionCmp,PeakDTOrganize, SIMPLIFY = F))))
    
  }
  
  PeakAnnotations<-mapply(function(x,y) {as.list(dt)[[3]][x:y]},
                          x= c(0, cumsum(NumPeaks[-length(Names)]))+1 ,y= cumsum(NumPeaks), SIMPLIFY = F)
  
  PeakAnnotations<-  lapply(PeakAnnotations, function(x) { paste(x, collapse = ";")})
  
  
  if(massOffset==0){
    Tags<-paste0("mods:",ModString," ", "ions:", PeakAnnotations)
  } else {
    Tags<-paste0("MassOffset: ",massOffset," mods:",ModString," ", "ions:", PeakAnnotations )
  }
  # Tags<-(massOffset)            
  
  # Tags<-paste0("mods:",ModString," ", "ions:", PeakAnnotations )
  # } else {
  #   Tags<-paste0("MassOffset: ",massOffset,"mods:",ModString," ", "ions:", PeakAnnotations )
  # }
  parallelTable<- data.table(blobMass=blobMass, blobInt=blobInt, 
                             PrecursorMasses=PrecursorMasses,Names=Names,Tags=Tags, FragmentationMode=FragmentationMode,
                             CollisionEnergy=CollisionEnergy)
  
  return(parallelTable)
  
}

LibraryParserSpxCmp<-cmpfun(SpXLibraryParser)

