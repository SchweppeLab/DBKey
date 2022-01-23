library(Rcpp)
sourceCpp('string_split.cpp')

testVector <- rep('I~am~a~boy',10)
for (i in 1:4){
  print(string_split(testVector,'~',1,i))
}



LibraryParser <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutoff,massOffset) {
  library(stringr)                                                                                                                                                                                          
  library(stringi)                                                                                                                                                                                          
  library(data.table,warn.conflicts = FALSE )                                                                                                                                                                                       
  suppressPackageStartupMessages(library(purrr))                                                                                                                                                                                            
  library(blob)
  library(readr)
  
  
  PrositLib<- (Library)
  
  nameindexes<-c(which(stri_detect_fixed(PrositLib,"Name: ")))
  headerLength<- which(stri_detect_fixed(PrositLib[1:100],"peaks:"))[1]-nameindexes[1]
  peakindexes<-nameindexes+headerLength
  HeaderLists<- unlist(mapply(function(x, y) {PrositLib[x:y]}, x = nameindexes, y = peakindexes, SIMPLIFY = T))
  PeakLists<- mapply(function(x, y) {PrositLib[x:y]}, x = peakindexes+1, y = c(nameindexes[-1]-1, length(PrositLib)))
  
  
  
  strt<-Sys.time()
  # PrositLib<-stri_replace_all_fixed(PrositLib,"FullName: ", "")
  # PrositLib<-stri_replace_all_fixed(PrositLib,"AvePrecursorMz:", "")
  PrositLib<-gsub("FullName: ", "", HeaderLists, fixed = T)
  PrositLib<-gsub("AvePrecursorMz:", "", HeaderLists, fixed = T)
  
  #getting library entry names, should be in "PEPTIDEK/2" format
 Names<- HeaderLists[stri_detect_fixed(HeaderLists, "Name: ")]
 # Names <- which(grepl("Name: :", PrositLib, fixed = T))
  Names <- str_remove_all(Names, "Name: ")
  #end<-Sys.time()
  #difftime(end,strt)
  #Names without charges
  #strt<-Sys.time()
  
  PeptideSequence <- gsub('.{2}$', '', Names, perl = T)
  
  #Retrieve mod string for further processing 
  ModString<- str_extract(HeaderLists,"ModString=[^=]+")
  ModString<- ModString[!is.na(ModString)]
  ModString<- str_extract(ModString,"//[^=]+")
  ModString<- gsub('.{6}$', '', ModString, perl = T)
  ModString<- gsub('^.{2}', '', ModString, perl = T)
  ModString<- gsub(' ', '', ModString, perl = T)
  
#  end<-Sys.time()
  #difftime(end,strt)
  
  #testing until I figure out how to read unimod xml
  #unimodTable <- read.csv("~/Repos/MSPtoDB/testMods.csv")
  #ModString<-stri_replace_all_fixed(ModString, unimodTable$mod, as.character(unimodTable$massshift))
  
  
  #Getting charge state + last aa for adding mods
 # strt<-Sys.time()
  
  Charge<-as.numeric(str_sub(Names,-1 ))
  LastAA<-str_sub(Names,-3,-3 )
  isK <- if(Source == "Prosit"){
    LastAA=="K"
  } else {
    LastAA=="]"
  }
  

  #Getting list of precursor masses, listed differently between Prosit/SpectraST
  PrecursorMasses <- if(Source == "Prosit"){
    PrositLib[stri_detect_fixed(HeaderLists, "MW:")]
  } else {
    PrositLib[stri_detect_fixed(HeaderLists, "PrecursorMZ:")]
  }
  PrecursorMasses <- gsub("[[:alpha:]]|:", "", PrecursorMasses)
  
  
  #If adding TMTPro, add mass to precursor mass
  PrecursorMasses<- if(TMTPro == TRUE) {
    as.numeric(PrecursorMasses)+304.207146/Charge
  } else  {PrecursorMasses}
  
  PrecursorMasses[isK] <- if(TMTPro == TRUE) {
    PrecursorMasses[isK]+304.207146/Charge[isK]
  } else {  PrecursorMasses[isK]}
  end<-Sys.time()
  difftime(end,strt)
  
  strt<-Sys.time()
  
  #Getting indicies of where peak lists are for further processing
  # peakindexes<-if(Source=="Prosit") {
  #   which(stri_detect_fixed(PrositLib,"peaks:"))+1
  # } else {
  #   which(stri_detect_fixed(PrositLib,"Peaks:"))+1
  # }
PeakLists[1]
  
   #I cant remember what problem I was trying to solve with this so I dont want to delete it yet
#  nameindexes<-c(which(stri_detect_fixed(PrositLib,"Name: "))[-1]-1,length(PrositLib))
  #cut<-min(length(peakindexes),length(nameindexes),length(PrecursorMasses))
  # nameindexes<-nameindexes[1:cut]
  # peakindexes<-peakindexes[1:cut]
  # PrecursorMasses<-PrecursorMasses[1:cut]
  # Names<-Names[1:cut]
  # LastAA<-LastAA[1:cut]
  # isK<-isK[1:cut]
  end<-Sys.time()
  difftime(end,strt)
  
  # 
  #retrieve peak lists from the indicies
  strt<-Sys.time()
  
  # PeakLists<- (mapply(function(x, y) {PrositLib[x:y]}, x = peakindexes, y = nameindexes))
  # PeakLists3<- (Map(function(x, y) {PrositLib[x:y]}, x = c(1,nameindexes[-length(nameindexes)]+1), y = peakindexes-1))
  
  PeakLists1<-mapply(function(x) {str_split(PeakLists[[x]], "\\t",simplify = T)}, x = seq(from=1,to=length(PeakLists), by=1))
  #PeakLists1<-mapply(function(x) {str_split(PeakLists[[x]], "\\t",simplify = T)}, x = seq(from=1,to=length(PeakLists), by=1))
  
   #PeakLists2<-mapply(function(x) {string_split(PeakLists[[x]], "\t",1,0)}, x = seq(from=1,to=length(PeakLists), by=1))
  
  #names(PeakLists)<-PeptideSequence 
  end<-Sys.time()
  difftime(end,strt)
  
  
  
  bYFunction<-function(x){
    pep<-unlist(str_split(x,""))
    b<- seq(1,length(pep), by=1)
    return(data.frame(aa=pep, b=paste0("b",b),y=paste0("y",y)))
  }
  
  
  # OrganizePeaks<- function(x) {
  #  dt<-data.table(masses=x)
  #  dt<- dt[, c("masses", "int", "annotaitons") :=(tstrsplit(masses, "\t"))]
  #  dt$masses<- dt[,substring(masses,2)]
  #  dt$masses<-as.numeric(dt$masses)
  #  dt$int<-as.numeric(dt$int)
  # 
  #   dt<-dt[dt$int>0,]
  #  # maxPeak<-max(dt$int)
  #   dt<-dt[(dt$int/maxPeak)>cutoff,]
  #   if(topX < length(dt$masses)) {
  #     peakNum <- min(topX, length(dt$masses))
  #      dt<-setorder(dt, -int)
  #     dt <-dt[1:peakNum,]
  #   }
  #   #dt<-setorder(dt, masses)
  #   return(dt)
  # }
  
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
  
  strt<-Sys.time()
  
  PeaksDT<- lapply(PeakLists, function(x){OrganizePeaks(x)}) 
  #OrganizePeaks(PeakLists[[1]])
  end<-Sys.time()
  difftime(end,strt)
  
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
  #Library<<-NULL
  return(parallelTable)
  
}
