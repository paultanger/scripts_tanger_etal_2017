# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

load("DS2013AllMeansNewWithCTD.Robject")

DS2013AllMeans$phenotype = as.factor(DS2013AllMeans$phenotype)
levels(DS2013AllMeans$phenotype)

# change order
DS2013AllMeans$phenotype  <- factor(DS2013AllMeans$phenotype , levels = 
                                      c("DTF",
                                        "DTH",
                                        "biomass",
                                        "height",
                                        "grain",
                                        "HI",
                                        "CTD",
                                        "Chla",
                                        "HTPheight",
                                        "NDRE",
                                        "NDVI",
                                        "NIR_760",
                                        "R_670",
                                        "RE_730"))

# define final phenos to include with all
finalphenos = c("DTF", "biomass", "height", "grain", "HI")
finalphenos = subset(DS2013AllMeans, phenotype %in% c(finalphenos))
finalphenos = ""
CanTemp   = subset(DS2013AllMeans, phenotype %in% c("CTD", finalphenos))
Chla      = subset(DS2013AllMeans, phenotype %in% c("Chla", finalphenos))
HTPheight = subset(DS2013AllMeans, phenotype %in% c("HTPheight", finalphenos))
NDRE      = subset(DS2013AllMeans, phenotype %in% c("NDRE", finalphenos))
NDVI      = subset(DS2013AllMeans, phenotype %in% c("NDVI", finalphenos))
NIR_760   = subset(DS2013AllMeans, phenotype %in% c("NIR_760", finalphenos))
R_670     = subset(DS2013AllMeans, phenotype %in% c("R_670", finalphenos))
RE_730    = subset(DS2013AllMeans, phenotype %in% c("RE_730", finalphenos))

# make a list
phenotypes = list(CanTemp, Chla, HTPheight, NDRE, NDVI, NIR_760, R_670, RE_730)

# name items in list
phenoNames = c("CTD", "Chla", "HTPheight", "NDRE", "NDVI", "NIR_760", "R_670", "RE_730")
names(phenotypes) = phenoNames

# drop levels
phenotypes = lapply(phenotypes, function(x) droplevels(x))

# for each item in list, make a new list with just the earliest useful data for that pheno
setwd(sourcepath)
source("GetEarlyPhenosFunction.R")
setwd(datapath)
earlyphenos = phenotypes
earlyphenos_rep1 = getearlydata1(earlyphenos)
earlyphenos_rep2 = getearlydata2(earlyphenos)
earlyphenos_rep3 = getearlydata3(earlyphenos)

# list of lists
earlyphenos = list(earlyphenos_rep1, earlyphenos_rep2, earlyphenos_rep3)
# name them
names(earlyphenos) = seq(1:3)
# save it
setwd(datapath)
filename = addStampToFilename("earlyphenosCTD", "Robject")
save(earlyphenos, file=filename)

# do the same for the overall best date
setwd(sourcepath)
source("GetBestPhenosFunction.R")
setwd(datapath)
bestphenos = phenotypes
bestphenos_rep1 = getBestdata1(bestphenos)
# bestphenos_rep2 = getBestdata2(bestphenos)
# bestphenos_rep3 = getBestdata3(bestphenos)


# save it
setwd(datapath)
filename = addStampToFilename("bestphenosCTD", "Robject")
save(bestphenos_rep1, file=filename)

# combine adjacent dates (best correlated) for each rep then get lsmeans

# combine sets of dates for 3 reps then get lsmeans

# keep each list, and also keep the final phenos list
filename = addStampToFilename("finalphenosDS2013", "Robject")
save(finalphenos, file=filename)

# what if we just keep ALL the phenotypes?
library(reshape2)
DS2013AllMeans$pheno_rep_date = paste(DS2013AllMeans$phenotype, DS2013AllMeans$rep, DS2013AllMeans$date, sep="_")

# rename line to ID
colnames(DS2013AllMeans)[1] = "ID"

# remove parents (not sure how they got there..)
parents = c("IR64", "Aswina")
DS2013AllMeans = DS2013AllMeans[!(DS2013AllMeans$ID %in% parents),]

# what should we keep
as.data.frame(colnames(DS2013AllMeans))

# make it wide
allphenos = dcast(DS2013AllMeans[,c(1,5,8)], ID ~ pheno_rep_date , value.var="mean", mean, na.rm=T)
# this might be a crazy idea.. 258 phenotypes?!
# but this will be way easier to compare things interactively in JMP with all QTL in one place

# save this for Rqtl in a place we can access via the server
setwd(datapath)
filename = addStampToFilename("AllDS2013HTPphenosByRepForRqtl", "csv")
write.table(allphenos, file=filename, na="-", col.names=T, row.names=F, sep=",", quote=FALSE)

# make qtl pipeline
