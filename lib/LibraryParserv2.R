
library(stringr)                                  
library(stringi)
library(data.table,warn.conflicts = FALSE )  
suppressPackageStartupMessages(library(purrr))  
library(blob)
library(readr)
library(compiler)


## Filtering peak functions
OrganizePeaks<- function(x,topX,cutoff,IonTypes) {
  
  dt<-x[which(x$int >0 ),]
  if(!is.null(IonTypes)) {
    dt<-dt[Reduce(`|`, lapply(IonTypes, grepl, x = dt$annotations)), ]
  }
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
  
  dt$annotations <- gsub("\"", "", dt$annotations)

  dt<-setkey(dt, masses)
  return(dt)
}


OrderPeaks<- function(x) {
  
  dt<-setkey(x, masses)
  return(dt)
}



OrgPeakCmp<-cmpfun(OrganizePeaks)
OrdPeakCmp<-cmpfun(OrderPeaks)

## Rebuilding function from OrgMassSpecR for monoisotopic mass correction https://orgmassspec.github.io/
ConvertPeptide <- function(sequence, output = "elements", IAA = FALSE) {

  peptideVector <- strsplit(sequence, split = "")[[1]]
  
  if(output == "elements") {
    
    FindElement <- function(residue) {
      
      if(residue == "A") element <- c(C = 3, H = 5, N = 1, O = 1, S = 0)
      if(residue == "R") element <- c(C = 6, H = 12, N = 4, O = 1, S = 0)
      if(residue == "N") element <- c(C = 4, H = 6, N = 2, O = 2, S = 0)
      if(residue == "D") element <- c(C = 4, H = 5, N = 1, O = 3, S = 0)
      if(residue == "E") element <- c(C = 5, H = 7, N = 1, O = 3, S = 0)
      if(residue == "Q") element <- c(C = 5, H = 8, N = 2, O = 2, S = 0)
      if(residue == "G") element <- c(C = 2, H = 3, N = 1, O = 1, S = 0)
      if(residue == "H") element <- c(C = 6, H = 7, N = 3, O = 1, S = 0)
      if(residue == "I") element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
      if(residue == "L") element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
      if(residue == "K") element <- c(C = 6, H = 12, N = 2, O = 1, S = 0)
      if(residue == "M") element <- c(C = 5, H = 9, N = 1, O = 1, S = 1)
      if(residue == "F") element <- c(C = 9, H = 9, N = 1, O = 1, S = 0)
      if(residue == "P") element <- c(C = 5, H = 7, N = 1, O = 1, S = 0)
      if(residue == "S") element <- c(C = 3, H = 5, N = 1, O = 2, S = 0)
      if(residue == "T") element <- c(C = 4, H = 7, N = 1, O = 2, S = 0)
      if(residue == "W") element <- c(C = 11, H = 10, N = 2, O = 1, S = 0)
      if(residue == "Y") element <- c(C = 9, H = 9, N = 1, O = 2, S = 0)
      if(residue == "V") element <- c(C = 5, H = 9, N = 1, O = 1, S = 0)
      
      if(residue == "C" && IAA == FALSE) element <- c(C = 3, H = 5, N = 1, O = 1, S = 1)
      if(residue == "C" && IAA == TRUE) element <- c(C = 5, H = 8, N = 2, O = 2, S = 1)
      
      return(element)
      
    }
    
    resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
    for(i in 1:length(peptideVector)) { resultsVector <- FindElement(peptideVector[i]) + resultsVector }
    
    resultsVector <- resultsVector + c(C = 0, H = 2, N = 0, O = 1, S = 0)   # add water
    
    return(as.list(resultsVector))
  }
  
  if(output == "3letter") {
    
    FindCode <- function(residue) {
      
      if(residue == "A") let <- "Ala"
      if(residue == "R") let <- "Arg"
      if(residue == "N") let <- "Asn"
      if(residue == "D") let <- "Asp"
      if(residue == "C") let <- "Cys"
      if(residue == "E") let <- "Glu"
      if(residue == "Q") let <- "Gln"
      if(residue == "G") let <- "Gly"	
      if(residue == "H") let <- "His"
      if(residue == "I") let <- "Ile"
      if(residue == "L") let <- "Leu"
      if(residue == "K") let <- "Lys"
      if(residue == "M") let <- "Met"
      if(residue == "F") let <- "Phe"
      if(residue == "P") let <- "Pro"
      if(residue == "S") let <- "Ser"
      if(residue == "T") let <- "Thr"
      if(residue == "W") let <- "Trp"
      if(residue == "Y") let <- "Tyr"
      if(residue == "V") let <- "Val"
      
      return(let)
      
    } 
    
    codes <- sapply(peptideVector, FindCode)
    return(paste(codes, collapse = ""))
    
  }
  
}

MonoisotopicMass <- function(formula = list(), isotopes = list(), charge = 0) {
  
  defaultFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0, Br = 0, Cl = 0, F = 0, Si = 0, x = 0)
  defaultFormula[names(formula)] <- formula   # replace default values with argument values
  
  defaultIsotopes <- list(C = 12, 
                          H = 1.0078250321, 
                          N = 14.0030740052, 
                          O = 15.9949146221, 
                          S = 31.97207069, 
                          P = 30.97376151,
                          Br = 78.9183376,
                          Cl = 34.96885271,
                          F = 18.99840320,
                          Si = 27.9769265327,
                          x = 0)
  
  defaultIsotopes[names(isotopes)] <- isotopes
  
  if(charge < 0 & abs(charge) > defaultFormula$H)
    stop("the number of negative charges exceeds the number of hydrogens in the formula list")
  
  mass <- (defaultFormula$C * defaultIsotopes$C + 
             defaultFormula$H * defaultIsotopes$H +
             defaultFormula$N * defaultIsotopes$N + 
             defaultFormula$O * defaultIsotopes$O +
             defaultFormula$S * defaultIsotopes$S + 
             defaultFormula$P * defaultIsotopes$P +
             defaultFormula$Br * defaultIsotopes$Br +
             defaultFormula$Cl * defaultIsotopes$Cl +
             defaultFormula$F * defaultIsotopes$F +
             defaultFormula$Si * defaultIsotopes$Si +
             defaultFormula$x * defaultIsotopes$x)
  
  if(charge != 0) mass <- abs((mass + charge * 1.007276466) / charge)
  
  return(mass)
  
}

LibraryParser <- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, CompoundClassArg,
                          Filter=FALSE, TMTPro=TMTPro, Source, topX=0, cutoff=0,massOffset=NA, IonTypes=NA) {

  # firstName<-grep("Name", Library[1:100])[1]
  # Library<-Library[firstName:length(Library)]
  
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
    match <- stri_detect_fixed(x, "ModString")
    if(any(match))
    {
      return(x[match])
    }
    else
    {
      return(x[which(stri_detect_fixed(x, "mods"))])
    }
  }
  Mods<-sapply(AltComments,stringFinder)
  
  unimodTable <- read.csv("~/Repos/MSPtoDB/unimod_custom.csv")
  unimodTable$mod <- str_escape(unimodTable$mod) #Ensure that characters won't interfere with regex replacement
  unimodTable$mod <- paste0("^", unimodTable$mod, "$") #Add start and end boundaries to treat as regex for exact match
  modforprecursor<-list()
  modparser <- function(x) { #Function takes in a single full ModString
    remove_modstring<-str_split(x,"//",simplify = T)[,2] #Remove the peptide sequence and "ModString=" content
    remove_modstring<-str_split(remove_modstring,"/",simplify = T)[,1] #remove the terminal charge "/2"

    if(nchar(remove_modstring)>=1){
      out<-lapply(str_split(remove_modstring,";",simplify = T), stringr::str_trim)
      out<-str_split(out,"@")
      mods<-as.data.frame(do.call(rbind,out))
      mod<-mods[,1]
      mod<-stri_replace_all_regex(mod, unimodTable$mod, as.character(unimodTable$massshift), vectorize_all = F)
      mod<-trimws(mod)
      modforprecursor<<-append(modforprecursor, list(mod))
      split_positions<-str_split(mods[,2], "[[:alpha:]]",simplify = T)
      number_cols<-min(2,ncol(split_positions))
      pos<-split_positions[,number_cols]
      pos<-str_split(pos, "/", simplify = T)[,1]
      if(number_cols<2)
      {
        # Use regular expression to extract the first continuous integer from the input
        for(i in 1:length(pos))
        {
          int_string <- regmatches(pos[i], gregexpr("[[:digit:]]+", pos[i]))[[1]]
          int_number <- as.numeric(int_string)
          # Decrement the number
          int_number <- int_number - 1

          # Replace the original integer with the decremented number in the input string
          pos[i] <- sub(int_string, as.character(int_number), pos[i], fixed=TRUE)
        }
      }
      pos[pos<=0 & pos!=""] <- 0
      returnstring<-pos
      returnstring[returnstring!=""] <- paste0(mod[returnstring!=""],"@",returnstring[returnstring!=""],";", collapse = "")
      returnstring<-trimws(returnstring,c("right"),";")
      returnstring<-trimws(returnstring,c("right")," ")
    } else {
      returnstring<-""
      modforprecursor<<-append(modforprecursor,c(0.0))
    }
    return(returnstring[[1]])
  }

  Modsoutput<-sapply(Mods,modparser)
  Modsoutput<-str_replace(Modsoutput," ","")
  has_rt<-sum(sapply(AltComments, function(x) {  return(stri_detect_regex(x,"RetentionTime|iRT"))})) # Check to see if iRT or RetentionTime is present.
  
  if(has_rt>=1)
  {
    rtItems<-sapply(AltComments, function(x) {  return(x[which(stri_detect_regex(x,"RetentionTime|iRT"))])    })
    RetentionTime <- str_split(rtItems, "=", simplify = T)[,2]
  } else
  {
      RetentionTime <- ""
  }
    
  getFrag<- function(x){
    match<- unique(stri_extract_all_fixed(x, c("CID","HCD"), simplify = TRUE, omit_no_match = TRUE))
    if(length(match)==2)
    {
      return(as.character(match[1]))
    }
    else
    {
      stop("Unable to determine Fragmentation method from file, please specify")
    }
  }
  
  if(FragmentationMode== "Read From file")
  {
    FragmentationMode=getFrag(HeaderLists[1:100])
  }
  else
  {
    FragmentationMode=FragmentationMode
  }
    

  getCE<- function(x){
    
    CE<-x[stri_detect_fixed(x, "Collision")]
    CE<-gsub("_", "", CE)
    CE<- stri_extract_first_regex(CE,"(?:Collisionenergy=)[\\d]+[\\.]?[\\d]+(?: )")
    CE<- gsub("Collisionenergy=", "", CE)
    CE<-trimws(CE,c("right")," ")
    if(length(CE)==0) {
      stop("Unable to determine Collision Energies from file, please specify")
    } else{
      return(CE)
    }
  }
      
    
if(CollisionEnergy== "Read from file")
    {
  CollisionEnergy=getCE(HeaderLists)
    }
    else
    {
      CollisionEnergy=CollisionEnergy
    }

  
  
  Names<- HeaderLists[stri_detect_fixed(HeaderLists, "Name: ")]
  Names <- str_remove_all(Names, "Name: ")
  
  
  NumPeaks<-(c(nameindexes[-1], length(Library)+1) - peakindexes)-1
  Charge<-as.numeric(str_sub(Names,-1 ))
  sequence<- gsub('.{2}$', '',Names)
    totalprecursormod<- lapply(modforprecursor, function(x) {sum(as.numeric(x))})
  totalprecursormod[is.na(totalprecursormod)] <- 0
  PrecursorMasses<- mapply(function(x,y,z) { 
    MonoisotopicMass(ConvertPeptide(x), charge = y) + (z/y)}, 
    x = sequence, y = Charge, z=totalprecursormod)



rm(HeaderLists)
  
  PeakLists <- mapply(function(x) { str_split(PeakLists[[x]], "\\s+", simplify = TRUE)},x = seq(from = 1, to = length(PeakLists), by=1))

  PeakLists <- na.omit(do.call("rbind", PeakLists))

  dt<-data.table(PeakLists)
  colnames(dt) <- c("masses", "int", "annotations")
  dt$masses<-as.numeric(dt$masses)
  dt$int<-as.numeric(dt$int)

  rm(PeakLists)
  if(sum(grepl( "/", dt$annotations, fixed = TRUE))>0)
  {
    dt$annotations<- gsub( "[^//]+$","",dt$annotations, perl = TRUE )
    dt$annotations<- gsub('[^\\^|[:alnum:]]', "", dt$annotations, perl=TRUE)
  }

    PeakDT<-mapply(function(x,y) {dt[x:y]},
          x= c(0, cumsum(NumPeaks[-length(Names)]))+1 ,y= cumsum(NumPeaks), SIMPLIFY = FALSE)

    PeakDTOrganize<- lapply(PeakDT, function(x){OrgPeakCmp(x,topX,cutoff,IonTypes)})
    rm(PeakDT)


    blobMass<-lapply(PeakDTOrganize, function(x) {(packBits(numToBits(unlist(x[,1]))))})
    blobInt<-lapply(PeakDTOrganize, function(x) {(packBits(numToBits(unlist(x[,2]))))})

    PeakAnnotations<-lapply(PeakDTOrganize, function(x) {
      paste(x$annotations,collapse=";")
    })
    
    
   
   

  if(length(massOffset) != 0)
  {
   seqCharge <- data.frame(tstrsplit(Names, "/"))
   seqCharge$mZ <- PrecursorMasses
   names(seqCharge) <- c("Sequence", "charge", "mZ")
   joined <-dplyr::left_join(seqCharge, read.csv(massOffset$datapath),by = "Sequence")
   joined$massOffsetTag<- ""
   joined$massOffsetTag[!is.na(joined$massOffset)] <- paste0("massOffset:",joined$massOffset[!is.na(joined$massOffset)])
    Tags<-paste0(joined$massOffsetTag," mods:",Modsoutput," ", "ions:", PeakAnnotations )
  }
  else
  {
    Tags<-paste0("mods:", Modsoutput, " ", "ions:", PeakAnnotations)  
  }

  # Handle mapping of compound class:
  mappedCompoundClasses = ""
  if(length(CompoundClassArg) != 0)
  {
   seqCharge <- data.frame(tstrsplit(Names, "/"))
   names(seqCharge) <- c("Sequence", "charge")
   joined <-dplyr::left_join(seqCharge, read.csv(CompoundClassArg$datapath),by = "Sequence")
   mappedCompoundClasses<-joined$CompoundClass
  }

  
     parallelTable<- data.table(blobMass=blobMass, blobInt=blobInt, 
                                 Formula =Modsoutput,
                             PrecursorMasses=PrecursorMasses,Names=Names,Tags=Tags, FragmentationMode=FragmentationMode,
                             CollisionEnergy=CollisionEnergy,CompoundClass=mappedCompoundClasses,
                             iRT=RetentionTime, seq=sequence, z= Charge)

  return(parallelTable)
  
}

LibraryParserCmp<-cmpfun(LibraryParser)


