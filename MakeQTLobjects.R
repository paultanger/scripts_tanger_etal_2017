# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

load("crossFinalLSmeansOrderFixedReRunCalcGeno_20150624_1158.Robject")
load("FinalLSmeansstepwiseDFafterfixchr_20150624_1332.RObject")

# remove epistatic
modeldf2 = lsmeansmodeldf[!is.na(lsmeansmodeldf$chromosome),]
# remove DTF
modeldf3 = modeldf2[!(modeldf2$phenotype == "DTF_LSmeans"),]

# modeldf3 = lsmeansmodeldf

# subset biomass
biomass = subset(modeldf3, phenotype == "biomass_LSmeans")
biomassmodel = makeqtl(cross, chr=biomass$chromosome, pos=biomass$position, qtl.name=biomass$id, what="prob")

biomassmodel = refineqtl(cross, pheno.col = "biomass_LSmeans", biomassmodel, method="hk")

# save these for lovell to play with
filename = addStampToFilename("biomassQTLmodel", "Robject")
# save(biomassmodel, file=filename)

plotLodProfile(biomassmodel)

load("biomassQTLmodel_20150727_1345.Robject")

library(RCurl)
script <- getURL("https://raw.githubusercontent.com/jtlovell/eQTL_functions/master/lodProfileGGplot.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# makeLpPlot(cross=crossObject, models=listOfModels, formulae=listOfFormulae, covar=nullUnlessCovariateWasUsed, allchr=T/F dependingOnIfYouWantAllChrsPlotted, ... PassedOnToGeom_line)

# make list of models
biomasslist = list(biomassmodel)
names(biomasslist) = "biomass_LSmeans"

# get list of formula from stepwiseqtl model
load("FinalLSmeansstepwiseoutputMaxQTL10AfterFixChrOrder_20150624_1227.Robject")
biomassformula = formula(modelList[[1]])
listOfFormulae = list(biomassformula)
names(listOfFormulae) = "biomass_LSmeans"
library(plyr)
makeLpPlot(cross=cross, models=biomasslist, formulae=listOfFormulae, covar=NULL, allchr=T)
plotLodProfile(modelList[[1]])

# all together
formulaList<-lapply(modelList, function(x) attr(x, "formula"))
makeLpPlot(cross=cross, models=modelList, formulae=formulaList, covar=NULL, allchr=F)
