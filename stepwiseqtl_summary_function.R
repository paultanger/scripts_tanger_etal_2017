# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

stepwiseqtl.summary<-function(cross, phes, max.qtls,pens, covar, printout, ci.method, drop=1, prob=.9){
  #set up the environment
  if(class(cross)=="riself" | class(cross)=="bc"){
    stats_out<-data.frame(phenotype=character(),chromosome=numeric(),position=numeric(),
                          df=numeric(),type3SS=numeric(),LOD=numeric(),perc.var=numeric(),Fstat=numeric(),P.chi2=numeric(),P.F=numeric(),
                          effect.estimate=numeric(),effect.SE=numeric(),effect.t=numeric(),
                          lowCImarker=character(), hiCImarker=character(),
                          lowCIpos=character(), hiCIpos=character())
  }else{
    stats_out<-data.frame(phenotype=character(),chromosome=numeric(),position=numeric(),
                          df=numeric(),type3SS=numeric(),LOD=numeric(),perc.var=numeric(),Fstat=numeric(),P.chi2=numeric(),P.F=numeric(),
                          lowCImarker=character(), hiCImarker=character(),
                          lowCIpos=character(), hiCIpos=character())
  }
  statcols<-names(stats_out)
  if(printout==T){cat("conducting stepwise QTL model selection","\n")}
  statsout<-data.frame()
  par(mfrow=c(4,1),mar = c(2,4,2,1))
  count=0
  model.out<-list()
  for (i in phes){
    #run the stepwise model
    count<-count+1
    if(printout==T){cat("phenotype (", count, ":", length(phes),") -- ",i,"\n", sep="")}
    stepout<-stepwiseqtl(cross,pheno.col=i, method="hk",
                         penalties=pens[i,], 
                         max.qtl=max.qtls, 
                         covar=covar,
                         refine.locations=T,additive.only=F, keeplodprofile=T,keeptrace=T,verbose=F)
    #skip the rest if it is a null qtl model
    if(nqtl(stepout)==0){
      if(printout==T){cat(paste(" ******* phenotype",i,"has a null qtl model","\n"))}
      next
    }
    model.out[[i]]<-stepout
    #plot the lod profile to the window
    plotLodProfile(stepout,main=paste(i,"formula: ", formula(stepout)))
    #fit the qtl model
    fit<-fitqtl(cross,
                pheno.col=i,
                qtl=stepout,
                formula=formula(stepout),
                get.ests=T,dropone=T,covar=covar,
                method="hk")
    #determine the number of terms/qtl in the model
    nqtls<-nqtl(stepout)
    nterms<-sum(countqtlterms(formula(stepout), ignore.covar=F)[c(1,4)])
    ncovar<-length(covar)
    #extract information if the qtl model has multiple qtls/covariates
    if(nterms>1){
      ests_out<-summary(fit)$ests[2:(nterms+1),1:3]
      cis<-data.frame()
      #calculate confidence intervals for each qtl
      for (j in 1:nqtls){
        if(ci.method=="drop"){ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=drop)}
        if(ci.method=="bayes"){ciout<-bayesint(stepout,qtl.index=j, expandtomarkers=F, prob=prob)}
        lowmarker<-rownames(ciout)[1]
        highmarker<-rownames(ciout)[3]
        lowposition<-ciout[1,2]
        highposition<-ciout[3,2]
        cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
      }
      chrs<-stepout$chr
      pos<-stepout$pos
      #if covariates, add in NAs into the dataframe
      if(ncovar>0) {
        covar.ci<-(rep(NA,4))
        cis<-rbind(covar.ci, cis)
        chrs<-c(names(covar), stepout$chr)
        pos<- c(names(covar), stepout$pos)
      }
      #if epistasis, add NAs into the dataframe
      if(as.numeric(countqtlterms(formula(stepout))[4])>0) {
        epi.ci<-data.frame()
        for (j in 1:countqtlterms(formula(stepout))[4]){
          epi.ci.out<-(rep(NA,4))
          cis<-rbind(cis,epi.ci.out)
        }
        chrs<-c(chrs,rep(NA,countqtlterms(formula(stepout))[4]))
        pos<- c(pos,rep(NA,countqtlterms(formula(stepout))[4]))
      }
      
      if(class(cross)=="riself" | class(cross)=="bc"){
        ests_out<-summary(fit)$ests[2:(nqtls+1),1:3]
        stats<-data.frame(rep(i,nterms),chrs,pos,fit$result.drop[1:nterms,],ests_out,cis)
        colnames(stats)<-statcols
      }else{
        stats<-data.frame(rep(i,nterms),chrs,pos,fit$result.drop[1:nterms,],cis)
        colnames(stats)<-statcols
      }
      
    }else{
      cis<-data.frame()
      for (j in 1:nqtls){
        if(ci.method=="drop"){ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=drop)
        }else{
          if(ci.method=="bayes"){ciout<-bayesint(stepout,qtl.index=j, expandtomarkers=F, prob=prob)}
        }
        lowmarker<-rownames(ciout)[1]
        highmarker<-rownames(ciout)[3]
        lowposition<-ciout[1,2]
        highposition<-ciout[3,2]
        cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
      }
      if(class(cross)=="riself" | class(cross)=="bc"){
        ests_out<-summary(fit)$ests[2:(nqtls+1),1:3]
        stats<-data.frame(c(i,stepout$chr[1],stepout$pos[1],
                            as.numeric(fit$result.full[1,c(1,2,4,5,3,6,7)]),as.numeric(ests_out),cis[1,]))
      }else{
        stats<-data.frame(c(i,stepout$chr[1],stepout$pos[1],
                            as.numeric(fit$result.full[1,c(1,2,4,5,3,6,7)]),cis[1,]))
      }
      ests_out<-summary(fit)$ests[2:(nqtls+1),1:3]
      
      colnames(stats)<-statcols; rownames(stats)<-stepout$name
    }
    stats_out<-rbind(stats_out,stats)
  }
  write.csv(stats_out, file="stepwisesummary.stats.csv", row.names=F)
  return(model.out)
}