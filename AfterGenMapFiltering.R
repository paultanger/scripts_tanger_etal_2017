# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

# save the raw qtldata object
setwd(datapath)
filename = addStampToFilename("rawqtldataobject99sim", "Robject")
save(qtldata, file=filename)

setwd(datapath)
filename = addStampToFilename("Fullpop_missingGenosAfterRemoveDupsAfterRoughMap", "png")
png(filename, width=60, height=8.5, units="in", res=150)
plot.missing(qtldata)
dev.off()

# this is sort of pointless since we made up the map just to import to Rqtl
filename = addStampToFilename("RoughMapAfterRemoveDups", "pdf")
pdf(filename, width=11, height=8.5)
plot.map(qtldata)
dev.off()

# check for seg distortion..well I think we already filtered those out?
genotable = geno.table(qtldata)
# filter for chi square test p value less than .0000001 - book example
segdismarkers = genotable[genotable$P.value < .00000001,]
filename = addStampToFilename("SegDist", "tsv")
write.table(segdismarkers, filename, col.names=T, row.names=F, sep="\t", quote=FALSE)

# check for identical individuals
cmpgeno = comparegeno(qtldata)
# reset plot parameters..
par(mfrow=c(1,1))
filename = addStampToFilename("IdenticalLines", "pdf")
pdf(filename, width=11, height=8.5)
hist(cmpgeno, breaks=100, xlab="prop of identical genotypes")
rug(cmpgeno)
dev.off()

# the above plot may be hard to read, so pull out the individuals (as pairs)
# this lists pairs of individuals with the same calls at more than 90% of markers
identicalgenos = which(cmpgeno > 0.99, arr.ind=T)
# save this
filename = addStampToFilename("IdenticalGenosList", "tsv")
write.table(identicalgenos, file=filename, col.names=T, row.names=F, sep="\t", quote=FALSE)

# remove one of each of these pairs of individuals
todrop = identicalgenos[,1]
qtldataremoveidentical = subset(qtldata, ind=!(qtldata$pheno$ID %in% todrop))
summary(qtldataremoveidentical)
summary(qtldata)

# remove high crossover lines
numcrossx = countXO(qtldataremoveidentical)
filename = addStampToFilename("CrossoversAfterRoughMap", "pdf")
pdf(filename, width=11, height=8.5)
plot(numcrossx, ylab="Number of crossovers")
dev.off()
summary(numcrossx)
mean(numcrossx[1:600])
mean(numcrossx)

# try to remove lines with high crossovers...
highcross = as.data.frame(numcrossx[numcrossx > 100])
todrop = as.numeric(row.names(highcross))
qtldata2 <- subset(qtldataremoveidentical, ind=!(qtldataremoveidentical$pheno$ID %in% todrop))
numcrossx = countXO(qtldata2)
summary(numcrossx)
filename = addStampToFilename("crossoversAFTERremoving5lines", "pdf")
pdf(filename, width=11, height=8.5)
plot(numcrossx, ylab="Number of crossovers")
dev.off()

# can get QTL in regions with low coverage.. not real.. so should check regions of low marker coverage..
# blue is entropy, red is variance.. should be zero..
filename = addStampToFilename("missingdata", "pdf")
pdf(filename, width=15, height=8.5)
plot.info(qtldata2, col=c("blue", "red"))
dev.off()

# plot histogram of individuals with missing markers
filename = addStampToFilename("histofmissingmarkerlines", "pdf")
pdf(filename, width=11, height=8.5)
hist(nmissing(qtldata2, what="mar"), breaks=50, ylab="number of markers", xlab="number of lines typed at markers")
dev.off()

# calculate genotyping errors.. LOD score of double crossovers..
# this takes a while... like 6 hours
# use haldane because that is what joinmap is using
ptm <- proc.time()
qtldata2 = calc.errorlod(qtldata2, map.function="haldane")
proc.time() - ptm

# save this
filename = addStampToFilename("qtldata2_afterRoughMapafterRemoveHighCrossesAfterCalcErrorLOD", "Robject")
save(qtldata2, file=filename)

# load data
# load("qtldata2_afterRoughMapafterRemoveHighCrossesAfterCalcErrorLOD_20150323_1318.Robject")
setwd("~/Desktop/Dropbox/working folder/Rprojects/RILpop/data/")

# show markers with highest LOD (maybe genotyped wrong)
top = top.errorlod(qtldata2, cutoff=5) # higher than LOD of 5
filename = addStampToFilename("LODofDoubleXOR", "tsv")
# write.table(top, file=filename, col.names=T, row.names=F, sep="\t", quote=FALSE)

# plot and remove potentially bad markers
allchr = names(qtldata2$geno)
filename = addStampToFilename("CheckHighLODmarkers", "pdf")
pdf(filename, width=11, height=8.5)
for (chr in allchr) {
  if (length(top$id[top$chr==chr]) > 1){
  plot.geno(qtldata2, chr, as.character(top$id[top$chr==chr]), cutoff=5)
  }
}
dev.off()

todrop = top
head(todrop)
# make these bad markers NA at these positions..
qtldata3 = qtldata2
# check that this works
row(qtldata3$pheno)[qtldata3$pheno$ID==1092]
qtldata3$pheno$ID==459
qtldata2$geno[["3"]]$data[459, "S3_17565838"]

for(i in 1:nrow(todrop)) {
  chr <- todrop$chr[i]
  # since geno$data doesn't have the line #, need row number of the line
  id <- row(qtldata3$pheno)[qtldata3$pheno$ID==todrop$id[i]]
  mar <- todrop$marker[i]
  qtldata3$geno[[chr]]$data[id, mar] <- NA 
}
# check that it worked
summary(qtldata2)
summary(qtldata3)
qtldata2$geno[["12"]]$data[258, "S12_14827890"]
qtldata3$geno[["12"]]$data[258, "S12_14827890"]

# save this
filename = addStampToFilename("rawdata2objectformergingwithgendistb4droponemarker", "Robject")
save(rawdata2, file=filename)

# maybe it make sense to do this next step earlier (like before joinmap)
# but I'm doing it here for now..
# use the GenDistPkg to remove 99% similar markers
setwd("~/Desktop/gitstuff/geneticdistance/")
source("BuildPackageIfYouHaveR312.R")
library(GeneticDistance)
library(qtl)
# load data
load("qtldata3_afterremovingDoubleXOR_20150331_0746.Robject")

# get out genos
genos = pull.geno(qtldata3)

# change back to A and B
genos[genos == 1] = "A"
genos[genos == 2] = "B"
genos[is.na(genos)] = "-"

genos = as.data.frame(genos)

# add line names back
finalgenos = cbind(qtldata3$pheno[1], genos)

geno_joinmap = finalgenos

# transpose
geno.t = t(geno_joinmap)

# add col names
colnames(geno.t) <- geno_joinmap[,1]

# remove first row
geno.t = geno.t[-1,]

# make row.names into column
geno.t2 = as.data.frame(cbind(line = rownames(geno.t), geno.t)) 

# add dumb column for joinmap
geno_joinmap = cbind(marker = geno.t2[,1], class="(a,b)", geno.t2[,2:ncol(geno.t2)])

# sort by marker this doesn't work
# geno_joinmap = geno_joinmap[with(geno_joinmap, order(marker)), ]

# save this
prefix = paste("Fullpop_IR64B_ForRemoving99Sim_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
write.table(geno_joinmap, file=paste(prefix, "tsv", sep=".") , col.names=TRUE, row.names=F, sep="\t", quote=FALSE)

# load file
setwd("~/Desktop/newmachome/Desktop/Dropbox/working folder/Rprojects/RILpop/data/")
rawdata = readFile("Fullpop_IR64B_ForRemoving99Sim_20150331_0750.tsv")

# define calls
callist = c("A", "B", "", "-")

ptm <- proc.time()
duplicates = finddups(rawdata, callist, threshold=99)
proc.time() - ptm
# Remove duplicate markers from rawdata.
rawdata2 = removedups(rawdata, duplicates)

# save to import to Rqtl.. but we don't have gen distance..

# get individual names back
ind = as.data.frame(colnames(geno_joinmap)[3:ncol(geno_joinmap)])
rawdata3 = cbind(ind, rawdata2)
colnames(rawdata3)[1] = "ID"

# transpose it
rawdata3t = setNames(data.frame(t(rawdata3[,-1])), rawdata3[,1])

# make row.names into column
rawdata3t = as.data.frame(cbind(marker = rownames(rawdata3t), rawdata3t))

# get group and cM from qtldata3 
# do first element
map = as.data.frame(qtldata3$geno[[1]]$map)
map$group = names(qtldata3$geno[1])
# then add the rest
i = 1
for (i in 2:length(qtldata3$geno)) {
  tempdf = as.data.frame(qtldata3$geno[[i]]$map)
  tempdf$group = names(qtldata3$geno[i])
  map = rbind(as.matrix(map), as.matrix(tempdf)) 
}

# make a df and fix colnames
map = as.data.frame(map)
map$marker = row.names(map)
colnames(map)[1] = "cM"
map = map[,c(3,2,1)]

# merge based on marker
rawdata4 = merge(map, rawdata3t, by="marker")
# sort
rawdata5 = rawdata4[with(rawdata4, order(group, cM)), ]

# save and import into rqtl again
source('~/Desktop/Rgit/rilpopr/FormatForRQTL.R')
forRQTL = FormatForRqtl(rawdata5)

filename = addStampToFilename("ForDropOneMarker99Sim", "csv")
write.table(forRQTL, file=filename,  na="", col.names=TRUE, row.names=F, sep=",", quote=FALSE)

# check for bad markers using droponemarker code (if it changes group length)
# this takes a long time
# try with parallel.. each chr separately..
setwd("~/Rgit/rilpopr/")
source("ParallelDropOneMarker.R")

qtldata6 <- read.cross("csv", ".", "ForDropOneMarker99Sim_20150331_1355.csv", genotypes=c("A","B"), na.strings=c(NA, "-"), estimate.map=F)
class(qtldata)[1] <-"riself"
# fix the order of chromosomes
qtldata6$geno = qtldata6$geno[c(4:14,1:3)]
summary(qtldata6)

# save object
filename = addStampToFilename("qtldata6", "Robject")
save(qtldata6, file=filename)

# load data
load("qtldata6_20150331_1358.Robject")

# load lovells data to test (smaller dataset)
# library(qtl)
# setwd("~/Dropbox/Rprojects/RILpop/data/")
# load("lovells.Robject")
# allchr = names(lovells$geno)[1:2]
# data(hyper)
# allchr = names(hyper$geno)[1:2]
# summary(hyper)
# plot.map(hyper)

# setwd("~/Dropbox/Rprojects/RILpop/data/")
# load("qtldata3_afterremovingDoubleXOR_20150324_1243.Robject")
allchr = names(qtldata6$geno)

ptm <- proc.time()
# this returns a list of the worst markers (this only took 80 min after removing 99% sim)
qtldata6droponemarker = runDropOneMarker(qtldata6, allchr, "haldane")
proc.time() - ptm

setwd("~/Desktop/Dropbox/")
# this is a scanone object, and can't be saved like a text file..
filename = addStampToFilename("droponemarker99simresults", "Robject")
save(qtldata6droponemarker, file = filename)

load("droponemarker99simresults_20150331_1848.Robject")

filename = addStampToFilename("droponemarkerresultsremove99sim", "pdf")
pdf(filename, width=11, height=8.5)
plot(qtldata6droponemarker, lod=1, ylim=c(-100,0))
plot(qtldata6droponemarker, lod=2, ylab="Change in chr length (cM)")
dev.off()

# list worst
worstmarkers = summary(qtldata6droponemarker, lodcolumn=2)

for (chr in allchr) {
  print(paste(chr, round(max(qtldata6$geno[[chr]]$map)), sep="  "))
}

# drop worst
badmar <- rownames(summary(qtldata6droponemarker, lod.column=2))
qtldata7 <- drop.markers(qtldata6, badmar)
summary(qtldata6)
summary(qtldata7)
2231-2217

# save object
filename = addStampToFilename("qtldata7", "Robject")
save(qtldata7, file=filename)

# load object
load("qtldata7_20150331_1550.Robject")

forjoinmap = FormatForJoinMap(qtldata7)

individuals = ncol(forjoinmap)-2
individuals
loci = nrow(forjoinmap)
loci

# output for joinmap
filename = addStampToFilename("AfterAllFiltering", "tsv")
write.table(forjoinmap, file=filename , col.names=TRUE, row.names=F, sep="\t", quote=FALSE)
