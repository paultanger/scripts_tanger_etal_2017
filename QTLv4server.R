# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

finaldata <- read.cross("csvs", ".", "genotypematrixRestMapFinal_20150514_1239.csv", "DS2013FinalPhenosLSmeans_20150608_1536.csv", genotypes=c("AA","BB"), na.strings="-", estimate.map=F)
alldata <- read.cross("csvs", ".", "genotypematrixRestMapFinal_20150514_1239.csv", "AllDS2013HTPphenosByRepForRqtl_20150520_1811.csv", genotypes=c("AA","BB"), na.strings="-", estimate.map=F)

finaldata = alldata

class(finaldata)[1] <-"riself"
summary(finaldata)

cross = finaldata
# load data to make plots for biomass mtg
load("crossallphenosOrderFixedReRunCalcGeno_20150612_1559.Robject")

# try this to clean up genotypes (from lovell)
cross <- subset(cross, ind=names(ntyped(cross))[ntyped(cross)!=0])
summary(cross)

# check stuff
filename = addStampToFilename("missingdata", "png")
png(filename, width=60, height=8.5, units="in", res=150)
plot.missing(cross)
dev.off()

filename = addStampToFilename("geneticmap", "png")
png(filename, width=11, height=8.5, units="in", res=150)
plot.map(cross)
dev.off()

# check for identical individuals
cmpgeno = comparegeno(cross)
# reset plot parameters..
filename = addStampToFilename("identicalgenos", "pdf")
par(mfrow=c(1,1))
pdf(filename, width=11, height=8.5)
hist(cmpgeno, breaks=100, xlab="prop of identical genotypes")
rug(cmpgeno)
dev.off()

# plot histogram of individuals with missing markers
filename = addStampToFilename("histofmissingmarkerlines", "pdf")
pdf(filename, width=11, height=8.5)
hist(nmissing(cross, what="mar"), breaks=50, ylab="number of markers", xlab="number of lines typed at markers")
dev.off()

# get a list of phenotypes..
phenos = as.data.frame(colnames(cross$pheno), stringsAsFactors=F)
phenos
# don't include ID column..
phenos = phenos[-1,]

# quick scanone check
# calc prob at 1cM intervals
ptm <- proc.time()
cross = calc.genoprob(cross, map.function="kosambi", step=1)
proc.time() - ptm


method = "hk"
# method = "ehk"
# method = "imp"
nodes = 23
ptm <- proc.time()
scan1 = scanone(cross, pheno.col=phenos, method=method, n.cluster=nodes)
proc.time() - ptm

# export it
filename = addStampToFilename("scan1finalLSmeans", "tsv")
write.table(scan1, file=filename, row.names=F, sep="\t", quote=F)

# save object
filename = addStampToFilename("scan1finalLSmeans", "Robject")
save(scan1, file=filename)

# also save cross object
filename = addStampToFilename("crossallphenos", "Robject")
save(cross, file=filename)

ptm <- proc.time()
scan1imp = scanone(cross.imp, pheno.col=phenos, method=method, n.cluster=nodes)
proc.time() - ptm
scan1 = scan1imp

# get rid of sci not
options(scipen=999)


filename = addStampToFilename("scan1impFinalPhenos", "pdf")
pdf(filename, width=10, height=8)
i= 1
# need to skip chr and pos columns
for(i in i:(length(scan1)-2)){
#   print (i)
  # but lodcolumn starts with first pheno (col 3)
  plot(scan1, lodcolumn=i, col=c("blue"), lty=1,      
       main=paste0("scanone imp QTLs Harvest Phenos ", colnames(scan1)[i+2]), # so add 2 to get the right label
       ylab="LOD")
}
dev.off()

# for each pheno, plot all 3 reps together
as.data.frame(colnames(scan1))
scan1sort = scan1[,c(1:3,11,19,
                     4,12,20,
                     5,13,21,
                     6,14,22,
                     7,15,23,
                     8,16,24,
                     9,17,25,
                     10,18,26)]
# final phenos
scan1sort = scan1[,c(1:5,7:9,11:13)]

filename = addStampToFilename("scan1FinalPhenosRepTogether", "pdf")
pdf(filename, width=10, height=8)
i= 1
# need to skip chr and pos columns and count by 3s
for(i in seq(i, length(scan1sort)-2, by=3)){
#     print (i)
  # but lodcolumn starts with first pheno (col 3)
  plot(scan1sort, lodcolumn=c(i:(i+2)), col=c("blue", "red", "black"), lty=1:3,     
       main=paste0("Harvest scanone HK QTLs, ", strsplit(colnames(scan1sort)[i+2],"_")[[1]][1]), # so add 2 to get the right label
       ylab="LOD")
  legend("topright", legend=colnames(scan1sort)[(i+2):(i+4)], lty=1:3, lwd=2, col=c("blue", "red", "black"))
}
dev.off()

# grain and HI
rep3stuff = scan1[,c(1:2,10,6)]
filename = addStampToFilename("scan1FinalPhenosGrainHI", "pdf")
pdf(filename, width=10, height=8)
plot(rep3stuff, lodcolumn=c(1:2), col=c("blue", "red"), lty=1:2,     
     main=paste0("Harvest scanone HK QTLs, grain, HI"), # so add 2 to get the right label
     ylab="LOD")
legend("topright", legend=colnames(rep3stuff)[3:4], lty=1:2, lwd=2, col=c("blue", "red"))
dev.off()

# save object
filename = addStampToFilename("crossB4permsfinalLSmeans", "Robject")
save(cross, file=filename)

setwd(datapath)
# load("crossB4perms_20150514_1634.Robject")

load("crossallphenos_20150521_1444.Robject")
load("scan1all_20150520_1857.Robject")

# calculate permutations
library(parallel)
detectCores()
perms = 1000
method = "hk"
nodes = 23

# get a list of phenotypes..
phenos = as.data.frame(colnames(cross$pheno), stringsAsFactors=F)
phenos
# don't include ID column..
phenos = phenos[-1,]

ptm <- proc.time()
allphenospermsscan1 = scanone(cross, pheno.col=phenos, method=method, n.perm=perms, n.cluster=nodes, addcovar=NULL, verbose=T)
proc.time() - ptm

# save it
filename = addStampToFilename("finallsmeansphenosscan1perms", "Robject")
save(allphenospermsscan1, file=filename)

# get this error on server.. but works fine on macbook..
# <simpleError in unclass(x)[repl, lodcolumn, drop = FALSE]: (subscript) logical subscript too long>

# then get this:
#   Error in updateParallelRNG(n.cluster) : object '.Random.seed' not found
set.seed(4321)

# try again, back to same error..
# something with the R version, it works with R 3.1.2 and any version of Rqtl:
# https://groups.google.com/forum/#!topic/rqtl-disc/DTyaRFwq87I
# so load it..
# load("earlyqtl1kpermsscan1.RObject")

load("finallsmeansphenosscan1perms_20150608_1546.Robject")

# can plot specific phenotypes
earlyqtl1kpermsscan1[,1]
max(earlyqtl1kpermsscan1[,1])
# or just get max
summary(earlyqtl1kpermsscan1, alpha=0.001)

# save max perms to df
alpha = 0.05 # lovell and mckay say 0.05 is common
df = data.frame()
for (i in 1:length(colnames(allphenospermsscan1)))
{
  # put name in first col
  df[i,1] = colnames(allphenospermsscan1)[i]
  # put value in second col
  df[i,2] = max(allphenospermsscan1[,i])
  df[i,3] = summary(allphenospermsscan1[,i], alpha=alpha)[1]
}

# rename stuff
permsSummary = df
colnames(permsSummary) = c("phenotype", "max", paste0("alpha_", alpha))

# sort
permsSummary = permsSummary[with(permsSummary, order(phenotype)), ]

# save it
filename = addStampToFilename("finalphenosscan1permsSummaryalpha0.05", "csv")
write.table(permsSummary, file=filename, row.names=F, sep=",", quote=F)

# so now check for phenotypes with sig QTL to continue with scan2 and model selection
load("allphenosscan1perms_20150528_1507.Robject")

setwd(sourcepath)
source("s1.qtlcheck.function.R")
setwd(datapath)

# define allphes for function
allphes = as.data.frame(colnames(cross$pheno), stringsAsFactors=F)
allphes = allphes[-1,]
qtl.test<- s1.qtlcheck(s1=scan1, s1.perms=allphenospermsscan1, alpha=0.001)

# or do it based on raw lod max instead
qtlpheshighLOD = qtl.test[as.numeric(as.character(qtl.test$s1.raw.max)) > 10, ]

#check to see which phenotypes have QTL
# qtlphes<-as.character(qtl.test$phenotype[qtl.test$s1.raw.qtl==TRUE])

# then eliminate some more data:
# 3/4, 4/17, raw phenos
qtlpheshighLOD$pheno = sapply(strsplit(as.character(qtlpheshighLOD$phenotype), "_"), `[`, 1)
rawspectra = c("NIR", "R", "RE", "DTH")
qtlpheshighLODfinal = subset(qtlpheshighLOD, !(pheno %in% rawspectra))

split = strsplit(as.character(qtlpheshighLODfinal$phenotype), "_")
# keep the last element of each list
split2 = vector()
# this is not optimized and takes forever for some reason..
for (i in 1:length(split)){
  split2 = c(split2, split[[i]][length(split[[i]])])
}  

qtlpheshighLODfinal$date = as.factor(split2)
levels(qtlpheshighLODfinal$date)

baddates = c("041713", "030413")
qtlpheshighLODfinal = subset(qtlpheshighLODfinal, !(date %in% baddates))


filename = addStampToFilename("phenoswithscan1HighLODphenos", "Robject")
save(qtlpheshighLODfinal, file=filename)

load("phenoswithscan1QTLalpha0.1_20150528_1728.Robject") 

# run scan2 perms
perms.out<-list()
ptm <- proc.time()
for (i in phenos){
  print(i)
  perms.out[[i]]<-scantwo(cross, pheno.col=i, method=method, n.perm=perms, n.cluster=nodes, addcovar=NULL, verbose=T)
}
proc.time() - ptm

# save it
filename = addStampToFilename("finalphenosscan2perms", "Robject")
save(perms.out, file=filename)

# so this took this long using 23 cores on the server:
# user      system     elapsed 
# 12225492.24    33056.86   540166.46 
# or 6.25 days

load("crossFinalLSmeansOrderFixedReRunCalcGeno_20150624_1158.Robject")
summary(cross)
plot.pheno(cross, pheno.col=5)

h <- read.csv("GrainAndHIDS2013formanuscript_20150722_1133.csv", header=TRUE, row.names=1)
head(h)

crossH <- add.phenos(cross, h, index="ID")
summary(crossH)
plot.pheno(crossH, pheno.col=6) #grain
plot.pheno(crossH, pheno.col=7) #HI
plot.pheno(crossH, pheno.col=8) #group

library(parallel)
detectCores()
perms = 1000
method = "hk"
nodes = 23

# run scan2 perms
perms.out<-list()
ptm <- proc.time()
for (i in 6:7){
  print(i)
  perms.out[[i]]<-scantwo(crossH, pheno.col=i, method=method, n.perm=perms, n.cluster=nodes, addcovar=NULL, verbose=T)
}
proc.time() - ptm

# save it
filename = addStampToFilename("HI_yield_scan2perms", "Robject")
save(perms.out, file="HI_yield_scan2perms.Robject")

perms.out[6:7]

save(crossH, file="crossFinalLSmeansOrderFixedReRunCalcGeno_20151120_HIyield.Robject")
