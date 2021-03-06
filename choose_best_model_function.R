# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.


#function to select the best model from two summary.stepwiseqtl
# two inputs needed:
#1 the list of models with a NULL covariate
#2 the list of models with an additive covariate
choose.best.model<-function(
  model.out=model.out,
  model.out.cyto=model.out.cyto,
  nocovar.modelstats=nocovar.modelstats,
  cytocovar.modelstats=cytocovar.modelstats,
  output)
{
  plod.nocovar<-sapply(model.out, function(x) attr(x, "pLOD"))
  plod.cytocovar<-sapply(model.out.cyto, function(x) attr(x, "pLOD"))
  plod.nocovar.df<-as.data.frame(plod.nocovar)
  plod.cytocovar.df<-as.data.frame(plod.cytocovar)
  plod.out<-merge(plod.nocovar.df,plod.cytocovar.df, by=0, all=T)
  plod.out[is.na(plod.out)]<-0
  best.model.type<-ifelse(plod.out$plod.nocovar>plod.out$plod.cytocovar,"no.covar","add.covar")
  model.types<-data.frame(plod.out$plod.nocovar,plod.out$plod.cytocovar,best.model.type)
  rownames(model.types)<-plod.out$Row.names
  colnames(model.types)<-c("plod.nocovar", "plod.cytocovar", "best.model.type")
  no.models<-model.out[rownames(model.types)[model.types$best.model.type=="no.covar"]]
  cyto.models<-model.out.cyto[rownames(model.types)[model.types$best.model.type=="add.covar"]]
  bestmodels<-append(cyto.models,no.models)
  
  stats.withcovar<-cytocovar.modelstats[complete.cases(cytocovar.modelstats),] #change this name to whatever you named the output from summary.stepwiseqtl
  stats<-nocovar.modelstats[complete.cases(nocovar.modelstats),]#change this name to whatever you named the output from summary.stepwiseqtl
  cyto.stats<-stats.withcovar[stats.withcovar$phenotype %in% rownames(model.types)[model.types$best.model.type=="add.covar"],]
  nocov.stats<-stats[stats$phenotype %in% rownames(model.types)[model.types$best.model.type=="no.covar"],]
  cyto.stats$position<-as.character(cyto.stats$position)
  nocov.stats$position<-as.character(nocov.stats$position)
  beststats<-rbind(cyto.stats,nocov.stats)
  beststats<-data.frame(beststats)
  ifelse(output=="models",return(bestmodels),
         ifelse(output=="stats",return(beststats),return(model.types)))
}