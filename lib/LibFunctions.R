library(compiler)

OrganizePeaks<- function(x) {
  
  dt<-x[which(x$int >0 ),]
  
  if(topX < length(dt$masses)) {
    peakNum <- min(topX, length(dt$masses))
    # dt<-setorder(dt, -int)
    # dt <-dt[1:peakNum,]
    dt$rank<- frank(dt$int)
    dt<-dt[dt$rank<=topX,]
    dt[,rank:=NULL]
  }
  if(cutoff>0) {
    maxPeak<-max(dt$int)
    dt<-dt[(dt$int/maxPeak)*100>cutoff,]
  }
  if(!is.na(IonTypes)) {
    dt<-dt[stri_detect_fixed(x$annotations,IonTypes),]
  }
  dt<-setkey(dt, masses)
  return(dt)
}

OrderPeaks<- function(x) {
  
  dt<-setkey(dt, masses)
  return(dt)
}

blobMassFunction<- function(x){
  as_blob(packBits(numToBits(unlist(x[,1]))))}
blobIntFunction<- function(x){
  as_blob(packBits(numToBits(unlist(x[,2]))))}

blobMassSort<- function(x,y) {
  as_blob(packBits(numToBits(as.list(dt)[[1]][x:y])))}
blobIntSort<- function(x,y) {
  as_blob(packBits(numToBits(as.list(dt)[[1]][x:y])))}


OrgPeakCmp<-cmpfun(OrganizePeaks)
OrdPeakCmp<-cmpfun(OrderPeaks)
blobMassFunctionCmp<- cmpfun(blobMassFunction)
blobIntFunctionCmp<- cmpfun(blobIntFunction)
blobMassSortCmp<- cmpfun(blobMassSort)
blobIntSortCmp<- cmpfun(blobIntSort)

