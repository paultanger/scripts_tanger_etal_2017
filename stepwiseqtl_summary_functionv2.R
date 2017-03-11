# stepwiseqtl.summary
# R function
# requires- qtl
# Author: John T. Lovell
# Version: 3.2
# Date: 18-Sept 2014


#stepwiseqtl.summary function details
#inputs
# 1. cross- a R/qtl cross object w/ genotype probabilities (only for method="hk" right now)
# 2. phes- a vector of phenotype names that correspond to phenotype columns in the R/qtl cross object
# 3. max.qtls- the maximum number of QTLs to explore via forward model selection
# 4. pens- a dataframe of penalties with column names that match phes (from calc.penalties ideally- have had trouble with other formats)
# 5. covar- a dataframe of covariate(s) with a length that matches that of the phenotypes in the cross object.
#         - interactive covariates are not incorporated. I am working on a post-hoc script to do this
# 6. printout- should the function tell you what is going on?
# 7. drop- the LOD drop interval that should be created. See the R/qtl function lodint

#outputs:
# 1. the function returns a list of stepwiseQTL objects. This includes:
#   a) formula, b) qtl model, c) lodprofiles, d) stepwise trace
# 2. a text file (comma separated) called "stepwisesummary.stats.csv" that contains the following columns:
#   phenotpye, chromosome, position (cM), 
#   type III statistics for each qtl: df, sums of squares, LOD, %variance explained, F statistic, P-value (from chi2), P-value (from F)
#   effect t statistics: estimate, standard error, t. test statistic
#   confidence intervals by position and marker
#   ## Note: for covariates and epistasis, some of these columns return "NA"


#the function:
stepwiseqtl.summary<-function(cross, phes, max.qtls,pens, covar, printout, drop=1.5){
  #set up the environment
  stats_out<-data.frame(phenotype=character(),chromosome=numeric(),position=numeric(),
                        df=numeric(),type3SS=numeric(),LOD=numeric(),perc.var=numeric(),Fstat=numeric(),P.chi2=numeric(),P.F=numeric(),
                        effect.estimate=numeric(),effect.SE=numeric(),effect.t=numeric(),
                        lowCImarker=character(), hiCImarker=character(),
                        lowCIpos=character(), hiCIpos=character())
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
                get.ests=T,dropone=T,covar=covar)
    #determine the number of terms/qtl in the model
    nqtls<-nqtl(stepout)
    nterms<-length(summary(fit)$ests[,1])-1
    ncovar<-length(covar)
    #extract information if the qtl model has multiple qtls/covariates
    if(nterms>1){
      ests_out<-summary(fit)$ests[2:(nterms+1),1:3]
      cis<-data.frame()
      #calculate confidence intervals for each qtl
      for (j in 1:nqtls){
        ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=drop)
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
      #put together the multi-qtl stats dataframe
      if(as.numeric(countqtlterms(formula(stepout))[4])>0) {
        epi.ci<-data.frame()
        for (j in 1:countqtlterms(formula(stepout))[4]){
          epi.ci.out<-(rep(NA,4))
          cis<-rbind(cis,epi.ci.out)
        }
        chrs<-c(chrs,rep(NA,countqtlterms(formula(stepout))[4]))
        pos<- c(pos,rep(NA,countqtlterms(formula(stepout))[4]))
      }
      stats<-data.frame(rep(i,nterms),chrs,pos,fit$result.drop[1:nterms,],ests_out,cis)
      colnames(stats)<-statcols
    }
    #extract information if the qtl model has only a single qtl
    else{
      cis<-data.frame()
      for (j in 1:nqtls){
        ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=drop)
        lowmarker<-rownames(ciout)[1]
        highmarker<-rownames(ciout)[3]
        lowposition<-ciout[1,2]
        highposition<-ciout[3,2]
        cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
      }
      ests_out<-summary(fit)$ests[2:(nqtls+1),1:3]
      stats<-data.frame(c(i,stepout$chr[1],stepout$pos[1],
                          as.numeric(fit$result.full[1,c(1,2,4,5,3,6,7)]),as.numeric(ests_out),cis[1,]))
      colnames(stats)<-statcols; rownames(stats)<-stepout$name
    }
    stats_out<-rbind(stats_out,stats)
  }
  write.csv(stats_out, file="stepwisesummary.stats.csv", row.names=F)
  return(model.out)
}