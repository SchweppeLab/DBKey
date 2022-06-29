
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
    dt<-dt[!(substr(dt$annotations,1,1) %in% IonTypes)]
    
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
                          Filter=FALSE, TMTPro=FALSE, Source, topX=0, cutoff=0,massOffset=NA, IonTypes=NA) {

  nameindexes<-c(which(stri_detect_fixed(Library,"Name: ")))
  headerLength<- which(stri_detect_regex(Library[1:100],"(?i)peaks:"))[1]-nameindexes[1]
  peakindexes<-nameindexes+headerLength
  HeaderLists<- unlist(mapply(function(x, y) {Library[x:y]}, x = nameindexes, y = peakindexes, SIMPLIFY = TRUE))
  PeakLists<- mapply(function(x, y) {Library[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(Library)))
  
  HeaderLists<-gsub("FullName: ", "", HeaderLists, fixed = TRUE)
  HeaderLists<-gsub("AvePrecursorMz:", "", HeaderLists, fixed = TRUE)
  Comments<-HeaderLists[which(stri_detect_fixed(HeaderLists,"Comment: "))]
  AltComments <- str_split(Comments, " (?=\\w+=)")

  stringFinder <- function(x)
  {
    return(x[which(stri_detect_fixed(x,"ModString"))])
  }

  Mods<-sapply(AltComments,stringFinder)
  
  unimodTable <- read.csv("~/Repos/MSPtoDB/unimod_custom.csv")
  unimodTable$mod <- paste("\\b",unimodTable$mod,"\\b", sep="") #Add word boundaries to treat as regex for exact match

  modparser <- function(x) { #Function takes in a single full ModString
    remove_modstring<-str_split(x,"//",simplify = T)[,2] #Remove the peptide sequence and "ModString=" content
    remove_modstring<-str_split(remove_modstring,"/",simplify = T)[,1] #remove the terminal charge "/2"

    out<-str_split(remove_modstring,"; ",simplify = T)
    out<-str_split(out,"@")
    mods<-as.data.frame(do.call(rbind,out))
    if(length(mods)>=1){
      mod<-mods[,1]
      mod<-stri_replace_all_fixed(mod, unimodTable$mod, as.character(unimodTable$massshift), vectorize_all = F)
      mod<-trimws(mod)
      pos<-str_split(mods[,2], "[[:alpha:]]",simplify = T)[,2]
      pos<-str_split(pos, "/", simplify = T)[,1]
      pos[pos<=0 & pos!=""] <- 0
      returnstring<-pos
      returnstring[returnstring!=""] <- paste0(mod[returnstring!=""],"@",returnstring[returnstring!=""],";", collapse = "")
      returnstring<-trimws(returnstring,c("right"),";")
      returnstring<-trimws(returnstring,c("right")," ")
    } else {
      returnstring<-""
    }
    return(returnstring[[1]])
  }

  Modsoutput<-sapply(Mods,modparser)


  rtItems<-sapply(AltComments,function(x) {  return(x[which(stri_detect_regex(x,"RetentionTime|iRT"))])    })


  RetentionTime <- str_split(rtItems, "=", simplify = T)[,2]
    
    getFrag<- function(x){
     match<- unique(stri_extract_all_fixed(x, c("CID","HCD"), simplify = TRUE, omit_no_match = TRUE))
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

  #PeptideSequence <- gsub('.{2}$', '', Names, perl = TRUE)
  

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
  
  #   PrecursorMasses<- if(TMTPro == TRUERUE) {
  #   as.numeric(PrecursorMasses)+304.207146/Charge
  # } else  {PrecursorMasses}
  # 
  # PrecursorMasses[isK] <- if(TMTPro == TRUERUE) {
  #   PrecursorMasses[isK]+304.207146/Charge[isK]
  # } else {  PrecursorMasses[isK]}
  
rm(HeaderLists)
  
  PeakLists1<- mapply(function(x, y) {Library[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(Library)),SIMPLIFY = TRUE)
  PeakLists<-mapply(function(x) {str_split(PeakLists[[x]], "\\t",simplify = TRUE)}, x = seq(from=1,to=length(PeakLists), by=1))
  PeakLists <- do.call("rbind", PeakLists)

  dt<-data.table(PeakLists)
  colnames(dt) <- c("masses", "int", "annotations")
  dt$masses<-as.numeric(dt$masses)
  dt$int<-as.numeric(dt$int)
  
   UnSorted<- any(unlist(mapply(function(x,y) {is.unsorted(dt[[1]][x:y])}, 
          x= c(0, cumsum(NumPeaks[1:4]))+1 ,y= cumsum(NumPeaks[1:5]), SIMPLIFY = FALSE)))
  
  rm(PeakLists)
  
  dt$annotations<- gsub( "[^//]+$","",dt$annotations, perl = TRUE )
  dt$annotations<- gsub('[^\\^|[:alnum:]]', "", dt$annotations, perl=TRUE)
  

    {
    PeakDT<-mapply(function(x,y) {dt[x:y]},
          x= c(0, cumsum(NumPeaks[-length(Names)]))+1 ,y= cumsum(NumPeaks), SIMPLIFY = FALSE)

    PeakDTOrganize<- lapply(PeakDT, function(x){OrgPeakCmp(x,topX,cutoff,IonTypes)})
    rm(PeakDT)

    blobMass<-(as_blob(flatten(mapply(blobMassFunctionCmp,PeakDTOrganize, SIMPLIFY = FALSE))))

    blobInt<-(as_blob(flatten(mapply(blobIntFunctionCmp,PeakDTOrganize, SIMPLIFY = FALSE))))
    
    PeakAnnotations<-lapply(PeakDTOrganize, function(x) {
      paste(x$annotations,collapse=";")
    })
    
    
  } #else {
  #   PeakDT<-mapply(function(x,y) {dt[x:y]},
  #                  x= c(0, cumsum(NumPeaks[-length(Names)]))+1 ,y= cumsum(NumPeaks), SIMPLIFY = FALSE)
  # 
  #   PeakDTOrganize<- lapply(PeakDT, function(x){OrdPeakCmp(x)})
  #   rm(PeakDT)
  # 
  #   blobMass<-(as_blob(flatten(mapply(blobMassFunctionCmp,PeakDTOrganize, SIMPLIFY = FALSE))))
  # 
  #   blobInt<-(as_blob(flatten(mapply(blobIntFunctionCmp,PeakDTOrganize, SIMPLIFY = FALSE))))
  # 
  # }
  # 

  # if(length(massOffset$name > 1)){
  # seqCharge <- data.frame(tstrsplit(Names, "/"))
  # seqCharge$mZ <- PrecursorMasses
  # names(seqCharge) <- c("Sequence", "charge", "mZ")
  # joined <-dplyr::left_join(seqCharge, read.csv(massOffset$datapath),by = "Sequence")
  # joined$mZ[!is.na(joined$massOffset )] <- as.numeric(joined$mZ[!is.na(joined$massOffset )])+
  #   (joined$massOffset[!is.na(joined$massOffset )])/as.numeric(joined$charge[!is.na(joined$massOffset )])
  # 
  # PrecursorMasses <- joined$mZ
  # joined$massOffsetTag<- ""
  # joined$massOffsetTag[!is.na(joined$massOffset)] <- paste0("MassOffset: ",joined$massOffset[!is.na(joined$massOffset)])
  # Tags<-paste0(joined$massOffsetTag," mods:",ModString," ", "ions:", PeakAnnotations )
  # 
  # 
  # } else {
   # Tags<-paste0("mods:",ModString," ", "ions:", PeakAnnotations)
    Tags<-paste0("mods:", Modsoutput, " ", "ions:", PeakAnnotations)
  #}
     parallelTable<- data.table(blobMass=blobMass, blobInt=blobInt, 
                             PrecursorMasses=PrecursorMasses,Names=Names,Tags=Tags, FragmentationMode=FragmentationMode,
                             CollisionEnergy=CollisionEnergy,
                             RetentionTime=RetentionTime)

  return(parallelTable)
  
}

LibraryParserCmp<-cmpfun(LibraryParser)


