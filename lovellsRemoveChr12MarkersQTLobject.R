# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

load("qtldata.Robject")

# run est.rf and plot chr 11 and 12
qtldata = est.rf(qtldata)
plot.rf(qtldata, chr=c(11,12))

rfs<-pull.rf(qtldata, chr=c(11,12))

#slice through the bad markers on chr 12 from marker 40 on chr 11
slice1 <- rfs[grep("S12",rownames(rfs)),40]
plot(slice1) #cuttoff = .42

#see that the bad markers are below .42, put them into a list
badmars12<-names(slice1)[slice1<.42]
rfs2<-pull.rf(qtldata)

#slice through one of the offending chr 12 markers on the rest of the lgs, except 12
slice2<-rfs2[grep("S12",rownames(rfs2), invert=T),"S12_24457241"] 
plot(slice2) #slice is at 0.425

#see that the bad markers are below a threshold of 0.425
badmars.rest<-names(slice2)[slice2<.425]
#combined bad markers into a vector and drop them
allbadmars<-c(badmars.rest,badmars12)
qtldata.nobad<-drop.markers(qtldata, markers=allbadmars)
#see that the geno matrix w/o these markers is better
qtldata.nobad<-est.rf(qtldata.nobad)
plot.rf(qtldata.nobad)

sapply(chrnames(qtldata.nobad),function(x) plot.rf(qtldata.nobad, chr=x))

# save it
datapath = "~/Desktop/Dropbox/working folder/Rprojects/RILpop/data/"
sourcepath = "~/Desktop/Rgit/rilpopr/"

setwd(sourcepath)
source("myfunctions.R")
source("FormatForJoinMap.R")
setwd(datapath)
filename = addStampToFilename("qtldataAfterRemoveBadChr12", "Robject")
# save(qtldata.nobad, file=filename)

# extract markers for joinmap
forjoinmap = FormatForJoinMap(qtldata.nobad)
summary(qtldata.nobad)
individuals = ncol(forjoinmap)-2
individuals
loci = nrow(forjoinmap)
loci

# output for joinmap
filename = addStampToFilename("AfterRemovingChr12BadForJoinMap", "tsv")
write.table(forjoinmap, file=filename , col.names=TRUE, row.names=F, sep="\t", quote=FALSE)

