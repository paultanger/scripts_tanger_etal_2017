# john lovell
# Please cite appropriately.

setwd("../input")
#cullmarkers- a function to reduce the size of the cross object by retaining only the markers within a set
# cM window that has the least missing information
#this function takes a R/qtl cross object and returns a culled object with one marker per stepsize
#cross= the original cross object to cull
#stepsize= numeric- the size of windows to return a single marker

cullmarkers<-function(cross=cross, stepsize=1){
  par(mfrow=c(2,1))
  plot.map(cross, main="original map")
  best.mar.out<-vector()
  for(i in 1:length(chrnames(cross))){
    print(i)
    chr.in<-cross$geno[[i]]
    chr.in<-chr.in$data
    pos<-matrix(unlist(pull.map(cross)[i]))
    mar<-markernames(cross,chr=chrnames(cross)[i])
    chr.in.info<-data.frame(pos,mar)
    colnames(chr.in.info)[1]<-"pos"
    j<-0
    while(j<max(chr.in.info$pos)){
      mars<-chr.in.info[chr.in.info$pos>=j & chr.in.info$pos<j+stepsize,]
      j<-j+stepsize
      mar.names<-mars$mar
      if(length(mars[,1])>1){
        dat<-chr.in[,mars$mar]
        n.miss<-apply(dat,2,function(x) sum(is.na(x)))
        best.mar<-as.character(mar.names[n.miss==min(n.miss)])
        if(length(best.mar)>1){
          best.mar<-sample(best.mar,1)
        }
      }else{
        if(length(mars[,1])==1){
          best.mar<-as.character(mar.names)
        }else{next}
      }
      best.mar.out<-c(best.mar.out,best.mar)
    }
  }
  markers.to.drop<-markernames(cross)[!markernames(cross) %in% best.mar.out]
  cross2<-drop.markers(cross, markers.to.drop)
  plot.map(cross2, main="culled map")
  return(cross2)
}


#for example
#culledcross<-cullmarkers(cross=cross, stepsize=1)