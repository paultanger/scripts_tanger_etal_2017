# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

#function to check which phenotypes have qtl under which mapping conditions
#currently optimized for s1 w/ and w/o covariates and w/ and w/o transformed traits (called qn here)
s1.qtlcheck<-function(s1=s1,
                      s1.perms=s1.perms,
                      s1.cov=NULL,
                      s1.perms.cov=NULL,
                      s1.qn=NULL,
                      s1.perms.qn=NULL,
                      s1.qn.cov=NULL,
                      s1.perms.qn.cov=NULL,
                      alpha=0.1)
{
  s1.out<-data.frame()
  for (i in allphes){
    s1.thresh<-as.numeric(summary(s1.perms[,i], alpha=alpha))
    s1.max<-max(s1[,i])
    s1.qtl<-as.character(s1.max>=s1.thresh)
    s1.out<-rbind(s1.out,cbind(i,s1.thresh,s1.max,s1.qtl))
  }
  s1.test<-s1.out
  colnames(s1.test)<-c("phenotype","s1.raw.thresh","s1.raw.max","s1.raw.qtl")
  if(length(colnames(s1.cov)>0)){
    s1.out<-data.frame()
    for (i in allphes){
      s1.thresh<-as.numeric(summary(s1.perms.cov[,i], alpha=alpha))
      s1.max<-max(s1.cov[,i])
      s1.qtl<-as.character(s1.max>=s1.thresh)
      s1.out<-rbind(s1.out,cbind(s1.thresh,s1.max,s1.qtl))
    }
    colnames(s1.out)<-c("s1.raw.cov.thresh","s1.raw.cov.max","s1.raw.cov.qtl")
    s1.test<-cbind(s1.test,s1.out)
  }
  if(length(colnames(s1.qn)>0)){
    s1.out<-data.frame()
    for (i in allphes){
      s1.thresh<-as.numeric(summary(s1.perms.qn[,i], alpha=alpha))
      s1.max<-max(s1.qn[,i])
      s1.qtl<-as.character(s1.max>=s1.thresh)
      s1.out<-rbind(s1.out,cbind(s1.thresh,s1.max,s1.qtl))
    }
    colnames(s1.out)<-c("s1.qn.thresh","s1.qn.max","s1.qn.qtl")
    s1.test<-cbind(s1.test,s1.out)
  }
  if(length(colnames(s1.qn.cov)>0)){
    s1.out<-data.frame()
    for (i in allphes){
      s1.thresh<-as.numeric(summary(s1.perms.qn.cov[,i], alpha=alpha))
      s1.max<-max(s1.qn.cov[,i])
      s1.qtl<-as.character(s1.max>=s1.thresh)
      s1.out<-rbind(s1.out,cbind(s1.thresh,s1.max,s1.qtl))
    }
    colnames(s1.out)<-c("s1.qn.cov.thresh","s1.qn.cov.max","s1.qn.cov.qtl")
    s1.test<-cbind(s1.test,s1.out)
  }
  return(s1.test)
}
