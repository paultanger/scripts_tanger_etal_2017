
library("reshape2")
library("Hmisc")
library("lsmeans")
library("lme4")
library("plyr")
library("ggplot2")
library("grid")
library("gridExtra")
library("scales")
library("corrplot")
library("RColorBrewer")
library("qtl")
library("parallel")
library("doParallel")

# Define function for filenames with date stamp
addStampToFilename = function(filename, extension){
  filenameWithStamp = paste(filename, "_", format(Sys.time(),"%Y%m%d_%H%M"), ".", extension, sep="")
  return (filenameWithStamp)
}

# calculate AICc
getAICc = function(model){
  paramaters = attr(logLik(model), "df")
  observations = sapply(ranef(model)[1],nrow)
  AICc = AIC(model) + 2*paramaters*(paramaters+1)/(observations-paramaters-1) 
  return (AICc)
}

# return the response var from a model
getrespname = function(model){
  return (all.vars(terms(model
  ))[1])
}

# get variance components..
varcompsdf = function(model){
  varcomps = as.data.frame(print(VarCorr(model),comp=c("Std.Dev.","Variance")))
  varcomps = varcomps[,-c(2,3)]
  varcomps.component = varcomps[,1]
  varcomps.variance = varcomps[,2]
  varcomps.SD = varcomps[,3]
  varcompsdf = data.frame(varcomps.component, as.numeric(varcomps.variance), as.numeric(varcomps.SD))
  names(varcompsdf) = c("component", "variance", "SD")
  varcompsdf$varPCT = varcompsdf$variance/ sum(varcompsdf$variance)*100
  return(varcompsdf)
}