LibraryParser <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, TMTPro, Source, topX, cutoff,massOffset) {
  
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
