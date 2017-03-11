# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")
hap = read.delim(file="filtered.hmp.txt", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
colnames(hap)[colnames(hap)=="rs#"] = "marker"

# make sure they are still in chr & position order..
hap = hap[with(hap, order(chrom, pos)), ]

# remove some columns
hap <- hap[, c(1:4,12:ncol(hap))]

# make column names easier to read..
# get first element
split = strsplit(colnames(hap[5:ncol(hap)]), ":")
split = sapply(split, "[", 1)
# add beginning columns
split2 = c(colnames(hap[1:4]), split)
# use character vector as new col names
colnames(hap) = split2

######################################################################################

# convert hets to H
hap_with_H <- as.matrix(hap)

# try the long way... this takes forever..
# try converting to matrix..see above..
hap_with_H[hap_with_H == "M"] <- "H"
hap_with_H[hap_with_H == "R"] <- "H"
hap_with_H[hap_with_H == "W"] <- "H"
hap_with_H[hap_with_H == "S"] <- "H"
hap_with_H[hap_with_H == "Y"] <- "H"
hap_with_H[hap_with_H == "K"] <- "H"

# convert NA to blank..
hap_with_H[is.na(hap_with_H)] = ""
# and N to blank..
hap_with_H[hap_with_H == "N"] <- ""

######################################################################################
# get consensus for parental wells..

happarents = hap_with_H[,c(1,5:24)]

# transpose to count occurences..
happarentsT = t(happarents)
# replace col names
colnames(happarentsT) = happarentsT[1,]
happarentsT = happarentsT[-1,]

# define groups
groups = as.vector(c("oldasw", "oldasw", "oldasw",
                     "newasw", "newasw", "newasw", "newasw", "newasw", "newasw", "newasw",
                     "oldIR64", "oldIR64",
                     "newIR64","newIR64","newIR64","newIR64","newIR64","newIR64","newIR64","newIR64"))

happarentsT = cbind(groups, happarentsT)

# make it long?
# make sure you don't have reshape loaded!!
detach("package:reshape", unload=TRUE)
library(reshape2)
happarentsT = as.data.frame(happarentsT)
# warning ok:
#   attributes are not identical across measure variables; they will be dropped 
melted = melt(happarentsT, id.vars="groups")
counts = as.data.frame(ftable(melted))


# to wide
wide = dcast(counts, groups + variable ~ value, value.var = "Freq")
colnames(wide) = c("groups", "marker", "N", "A", "C", "G", "H", "T")

# sum of all calls
wide$sum = rowSums(wide[,c(3:8)])

# calculate percents
wide$Npct = wide$N/wide$sum
wide$Apct = wide$A/wide$sum
wide$Cpct = wide$C/wide$sum
wide$Gpct = wide$G/wide$sum
wide$Hpct = wide$H/wide$sum
wide$Tpct = wide$T/wide$sum

# keep the call with > 49%
# remove old cols
widefinal = wide[,c(1,2,10:15)]
colnames(widefinal) = c("group", "marker", "N", "A", "C", "G", "H", "T")

# make NA if < 0.5
widefinal[widefinal<0.5] = NA
# warning is ok

# make long and discard NA.. should have 1 call for each marker & group
# this probably also discards markers for which there are no good calls?
longfinal = melt(widefinal, id.vars=c("group","marker"), na.rm=T, value.name="freq")

# sort by group and marker
longfinal = longfinal[with(longfinal, order(group, marker)), ]

# check for duplicates.. there shouldn't be..unless there are 50/50 calls.. happens with new & old IR64.. since even numbers of wells..
# for newIR64, since we have 4 wells with a call, vs 4 wells with N, just discard the N rows..
longfinal2 = subset(longfinal, !(group == "newIR64" & variable == "N" & freq == 0.5))
# for oldIR64, harder to decide what to do..keep the one good call?
longfinal3 = subset(longfinal2, !(group == "oldIR64" & variable == "N" & freq == 0.5))

# make wide to merge with hap file
# test = dcast(longfinal3[,-4], marker ~ group, value.var="variable")
# there must be duplicates..
# so there are some ties between H and a good call in oldIR64...
# so discard these as well...
longfinal4 = subset(longfinal3, !(group == "oldIR64" & variable == "H" & freq == 0.5))

# ugh!! there are still some ties between different calls for oldIR64.. what a headache!
longfinal4$uniq = paste(longfinal4$group,longfinal4$marker, sep="_")
dups = longfinal4[duplicated(longfinal4$uniq),]

# so discard these too...
longfinal5 = longfinal4[!(duplicated(longfinal4$uniq) | duplicated(longfinal4$uniq, fromLast = TRUE)), ]

# make wide to merge with hap file
tomerge = dcast(longfinal5[,-4], marker ~ group, value.var="variable")

# so merge parents back into full list..
allmerged = merge(tomerge, hap_with_H, all=T)

# rearrange and sort
allmerged = allmerged[,c(1,6:11,4,12:18,2,19:20,5,21:28,3,29:ncol(allmerged))]

# just keep calls with good calls in both parents
tofilter = allmerged[,-c(5:15,17:27)]

# only keep if have a good call
goodcalls <- c("A", "C", "G", "T")
tofilter2 = tofilter[tofilter$newasw %in% goodcalls & tofilter$newIR64 %in% goodcalls,]
# only keep if polymorphic
tofilter3 = tofilter2[tofilter2$newasw != tofilter2$newIR64,]

hap_with_H = tofilter3

#######################################
#######################################
# filter markers first...
# so only keep markers with calls in 60% of population...

# this is how many lines with no call
countmiss = apply(hap_with_H[7:(ncol(hap_with_H))] == "", 1, sum)
# this is number of lines
allcalls = dim(hap_with_H)[2] - 6
# so divide lines with no call by total number of lines
# gives percent missing...
pctmissing = as.matrix(countmiss / allcalls)

# so keep markers with < 0.40 missing
hap_with_H = hap_with_H[pctmissing < 0.40, ]

######################################################################################
# filter by het freq...

# this calculates het frequency
het_freq <- as.matrix(apply(hap_with_H=="H",2,sum)/(nrow(hap_with_H)))

# filter large % hets
# this works by creating a new dataframe based on the old one 
#and evaluating if het_freq is true or false for each row.. and only keeping true
hapFINAL <- hap_with_H[,het_freq<0.15]


######################################################################################
### try to call based on IR64 as parent###

# first change het calls in parents to N
# this creates column arrays.. that line up with the matrix of calls
# we will use this to compare each column to..

colnames(hapFINAL)[colnames(hapFINAL)=="newasw"] = "Aswina"
colnames(hapFINAL)[colnames(hapFINAL)=="newIR64"] = "IR64"

Aswina <- hapFINAL$Aswina
IR64   <- hapFINAL$IR64

# now change H calls to N
Aswina[Aswina=="H"] = "N"
IR64[IR64=="H"] = "N"

# works MUCH faster as a matrix..
hapFINAL2 = as.matrix(hapFINAL)

# make a new matrix, same size as input file, to fill with A, B or - but fill with NA first
genotypes2 <- hapFINAL2
genotypes2[,c(7:ncol(genotypes2))] = NA

# loop through each row, checking if equal to IR64
# these loops work by looping over col as i... [row,col]
# and it updates the actual genotype file as it goes..
for(i in 7:ncol(hapFINAL2)){
  genotypes2[(hapFINAL2[,i] == IR64),i] = "B"  ## assuming IR64 is correct
  genotypes2[(hapFINAL2[,i] == Aswina),i] = "A"
  genotypes2[(hapFINAL2[,i]==""),i] = "-"
  genotypes2[(hapFINAL2[,i]=="H"),i] = "-"
}

fileprefix = paste("genotypesIR64asB_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(genotypes2, file=paste(fileprefix, "tsv", sep="."), sep="\t", quote=F, row.names=F, col.names=T) 

# remove parents
genotypes2 = genotypes2[,-c(5:6)]

# save as Robject
filename = addStampToFilename("genotypes2", "Robject")
# save(genotypes2, file=filename)

# back to dataframe
genotypes2 = as.data.frame(genotypes2)

# these are not sorted correctly..
genotypes2 = genotypes2[with(genotypes2, order(chrom, pos)), ]

# possibly here is where we want to remove duplicate markers?
# use the GenDistPkg to remove 99% similar markers
# this is a package from a manuscript in prep
setwd("geneticdistance/")
source("BuildPackageIfYouHaveR312.R")
library(GeneticDistance)

# need to format and save for package
rawdata = genotypes2[,-c(2:4)]

# save this
prefix = paste("Fullpop_IR64B_ForRemoving99Sim_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")

# save this
prefix = paste("Fullpop_IR64B_Genotypes2_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(rawdata, file=paste(prefix, "tsv", sep=".") , col.names=TRUE, row.names=F, sep="\t", quote=FALSE)

# load file
rawdata = readFile("Fullpop_IR64B_Genotypes2_20150330_1617.tsv", classification=F)

# define calls
callist = c("A", "B", "", "-")

ptm <- proc.time()
duplicates = finddups(rawdata, callist, threshold=99)
proc.time() - ptm
# Remove duplicate markers from rawdata.
rawdata2 = removedups(rawdata, duplicates)

# get individual names back
ind = as.data.frame(colnames(genotypes2)[5:ncol(genotypes2)])
rawdata3 = cbind(ind, rawdata2)
colnames(rawdata3)[1] = "ID"

# transpose it
rawdata3t = setNames(data.frame(t(rawdata3[,-1])), rawdata3[,1])

# make row.names into column
rawdata3t = as.data.frame(cbind(marker = rownames(rawdata3t), rawdata3t))

# get stuff from genotypes2 that we lost, match based on marker
rawdata4 = merge(genotypes2[,c(1:4)], rawdata3t, by="marker")
genotypes2 = rawdata4

# filter again for lines with lots of missing..first calculate it..
missing_freq <- as.matrix(apply(genotypes2=="-",2,sum)/(nrow(genotypes2)))

# filter them out.. with more than 20% missing
genotypes2 <- genotypes2[,missing_freq<0.2]
missing_freq <- as.matrix(apply(genotypes2=="-",2,sum)/(nrow(genotypes2)))
genotypes2 <- genotypes2[,missing_freq<0.15]

###########################
# save for later so we don't have to do this again..
###########################
setwd("~/Desktop/Dropbox/working folder/rprojects/RILpop/data/")
fileprefix = paste("genotypesfilteredlines_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(genotypes2, file=paste(fileprefix, "tsv", sep="."), sep="\t", quote=F, row.names=F, col.names=T) 

##########################

# filter for allele freq
# add up B calls, calculate % of total..
# later I should clean this up and keep the calculations in a separate df.... then I don't have to delete them or account for them..
genotypes2$Bsum = apply(genotypes2[5:(ncol(genotypes2))]=="B", 1, sum)
genotypes2$Asum = apply(genotypes2[5:(ncol(genotypes2)-1)]=="A", 1, sum)
genotypes2$notN = apply(genotypes2[5:(ncol(genotypes2)-2)]!="-", 1, sum)
genotypes2$Bfreq = genotypes2$Bsum / genotypes2$notN

# remove marker if not between 30%-70%
genotypes2 <- genotypes2[genotypes2$Bfreq>0.3 & genotypes2$Bfreq<0.7,]

# remove stuff
genotypes2[c(1:5),c(1000:1009)]
genotypes2 = genotypes2[,-c(ncol(genotypes2):(ncol(genotypes2)-3))]

# now that we have less lines, redo the marker filtering...

# this is how many lines with no call
countmiss = apply(genotypes2[5:(ncol(genotypes2))] == "-", 1, sum)
# this is number of lines
allcalls = dim(genotypes2)[2] - 4
# so divide lines with no call by total number of lines
# gives percent missing...
pctmissing = as.matrix(countmiss / allcalls)

# so keep markers with < 0.20 missing
genotypes3 = genotypes2[pctmissing < 0.20, ]
genotypes4 = genotypes3[!(duplicated(genotypes3[,5:ncol(genotypes3)])),]

# format for Rqtl with fake cM data
genotypes5 = formatforRqtlfakecM(genotypes4)

# save it for Rqtl
genoprefix = paste("PolandFullPopJoinmapMLdata_ForRqtl_B4newRfiltering_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(genotypes5, file=paste(genoprefix, "csv", sep="."),  na="", col.names=TRUE, row.names=F, sep=",", quote=FALSE)

genotypes5 = read.table("PolandFullPopJoinmapMLdata_ForRqtl_B4newRfiltering_20150330_1726.csv", sep=",", header=T)

# load in Rqtl
library(qtl)
qtldata <- read.cross("csv", ".", "PolandFullPopJoinmapMLdata_ForRqtl_B4newRfiltering_20150330_1726.csv", genotypes=c("A","B"), na.strings=c(NA, "-"), estimate.map=F)
# warning about big cM is ok.. fake anyways..
class(qtldata)[1] <-"riself"
summary(qtldata)

# check for identical individuals
cmpgeno = comparegeno(qtldata)

# pull out the individuals (as pairs)
# this lists pairs of individuals with the same calls at more than 90% of markers
identicalgenos = which(cmpgeno > 0.99, arr.ind=T)

# remove one of each of these pairs of individuals
todrop = identicalgenos[,1]
qtldataremoveidentical = subset(qtldata, ind=!(qtldata$pheno$ID %in% todrop))

# extract genotype data from object to use for joinmap, MSTmap, or our CS tool
genos = as.data.frame(qtldataremoveidentical$geno[[1]]$data)
i = 1
for (i in 2:length(qtldataremoveidentical$geno)) {
  genos = cbind(as.matrix(genos), qtldataremoveidentical$geno[[i]]$data) 
}

# just discovered that there is an easier way to do this:
# pull.geno()  and pull.pheno() also exists

# change back to A and B
genos[genos == 1] = "A"
genos[genos == 2] = "B"
genos[is.na(genos)] = "-"

genos = as.data.frame(genos)

# add line names back
finalgenos = cbind(qtldataremoveidentical$pheno[1], genos)

# save this
prefix = paste("Fullpop_IR64B_PreFiltering_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(finalgenos, file=paste(prefix, "tsv", sep=".") , col.names=TRUE, row.names=F, sep="\t", quote=FALSE)

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

# check stuff
individuals = ncol(geno_joinmap)-2
individuals
loci = nrow(geno_joinmap)
loci

# output for joinmap
prefix = paste("Fullpop_IR64B_ForJoinmap_PreFiltering_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(geno_joinmap, file=paste(prefix, "tsv", sep=".") , col.names=TRUE, row.names=F, sep="\t", quote=FALSE)