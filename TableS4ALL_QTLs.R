# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

load("modeldfHTP_20150616_1549.Robject")
modeldfHTP = modeldfHTP
# remove epistatic QTL
modeldfHTP = modeldfHTP[!is.na(modeldfHTP$chromosome),]
# change from factors to numbers
modeldfHTP$chromosome = as.character(modeldfHTP$chromosome)
modeldfHTP$lowCIpos = as.numeric(as.character(modeldfHTP$lowCIpos))
modeldfHTP$hiCIpos = as.numeric(as.character(modeldfHTP$hiCIpos))

# get final manual phenotype QTLs

# load("modeldfAllFinalPhenosWithGrainHI_20151127_1630.RObject")
# that is wrong
# try this.. it has grain and HI but the funky results for biomass, height, DTH
load("stepwiseDFGrainHI_20151127_1729.RObject")
# this also works I think but it doesn't have grain and HI..
load("FinalLSmeansstepwiseDFafterfixchr_20150624_1332.RObject")

# so extract grain and HI and add to the old file:
phes = c("grain_rep3", "HI_rep3")

modeldf = modeldf[(modeldf$phenotype %in% phes),]

# add to biomass, height, DTH
FINALmodeldf = rbind(lsmeansmodeldf,modeldf)

head(FINALmodeldf)
# remove epistatic
FINALmodeldf = FINALmodeldf[!is.na(FINALmodeldf$chromosome),]
# remove DTF
FINALmodeldf = FINALmodeldf[!(FINALmodeldf$phenotype == "DTF_LSmeans"),]

# merge them together
as.data.frame(colnames(modeldfHTP))
as.data.frame(colnames(FINALmodeldf))

# remove unnecessary columns
modeldfHTP = modeldfHTP[,-1]

# just keep HTP we want
htps<-c("NDRE", "NDVI",  "Chla", "HTPheight", "CTD")
modeldfHTP = modeldfHTP[modeldfHTP$phenotype %in% htps ,]

# do this in excel...
# export them

filename = addStampToFilename("HTP_QTL", "tsv")
write.table(modeldfHTP, file=filename, col.names=T, row.names=F, sep="\t", quote=FALSE)

filename = addStampToFilename("FinalQTL", "tsv")
write.table(FINALmodeldf, file=filename, col.names=T, row.names=F, sep="\t", quote=FALSE)

