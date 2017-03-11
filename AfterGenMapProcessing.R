# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

# in joinmap, create genetic maps using recomb freq in groupings tree, choose 14 groups, make maps using ML (defaults)
# on each map, export the list of markers with distance


# try to rbind all the files together
file_list = list.files()
rm(cMfromML)
cMfromML <- do.call("rbind",lapply(file_list,
                                  FUN=function(files){
                                    read.table(files, header=TRUE, sep="\t")}))

cMfromML$Group = as.factor(cMfromML$Group)
levels(cMfromML$Group)

genoprefix = paste("cMfromJoinmap_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(cMfromML, file=paste(genoprefix, "csv", sep="."),  col.names=TRUE, row.names=F, sep=",", quote=FALSE)

# load genotypes file to merge with cM
# use the input file from joinmap
setwd("~/Desktop/Dropbox/")
# finalgenos = read.table("Fullpop_IR64B_PreFiltering_20150320_1806.tsv", sep="\t", header=T)
joinmapfile = read.table("Fullpop_IR64B_ForJoinmap_PreFiltering_20150330_1731.tsv", sep="\t", header=T)

# merge together..
geno_rqtl_all = merge(cMfromML[,c(3,4,5)], joinmapfile, by.x="Locus", by.y="marker", all=T)
# change order of levels before sorting
geno_rqtl_all$Group = factor(geno_rqtl_all$Group, levels=c("1part1","1part2","2","3","4part1","4part2","5","6","7","8","9","10","11","12"))

# sort by group then position
geno_rqtl_all = geno_rqtl_all[with(geno_rqtl_all, order(Group, Position)), ]

# rename col
colnames(geno_rqtl_all)[colnames(geno_rqtl_all)=="Locus"] = "marker"
colnames(geno_rqtl_all)[colnames(geno_rqtl_all)=="Position"] = "cM"

# don't need classification column
geno_rqtl_all = geno_rqtl_all[,-4]

# rearrange
# geno_rqtl_all = geno_rqtl_all[,c(1,3,2,4:ncol(geno_rqtl_all))]
geno_rqtl = geno_rqtl_all

# transpose it
geno_rqtl.t = t(geno_rqtl)
head(geno_rqtl.t[1:5,1:5])

# add line col as "ID"
geno_rqtl.t <- cbind(ID = rownames(geno_rqtl.t), geno_rqtl.t)

# remove chrom and marker and cM from ID col
geno_rqtl.t[1,1] = ""
geno_rqtl.t[2,1] = ""
geno_rqtl.t[3,1] = ""
geno_rqtl.t[1:5,1:5]

# rename cols
colnames(geno_rqtl.t) <- as.character(geno_rqtl.t[1,])
# fix ID col
colnames(geno_rqtl.t)[1] = "ID"
# remove old row with col names
geno_rqtl.t = geno_rqtl.t[-1,]
str(geno_rqtl.t)

# remove X from IDs
split = strsplit(geno_rqtl.t[,1], "X")
split = sapply(split, "[", 2)
# use character vector as new col names
geno_rqtl.t = cbind(as.data.frame(as.character(split)), geno_rqtl.t)
geno_rqtl.t = geno_rqtl.t[,-2]
colnames(geno_rqtl.t)[colnames(geno_rqtl.t)=="as.character(split)"] = "ID"
geno_rqtl.t$ID = as.character(geno_rqtl.t$ID)
geno_rqtl.t[1,1] = ""
geno_rqtl.t[2,1] = ""

# try with numeric instead of character.. na as ""
geno_rqtl.t$ID = as.numeric(geno_rqtl.t$ID)
genoprefix = paste("AfterJoinMapRemoving99SimCleaned", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
write.table(geno_rqtl.t, file=paste(genoprefix, "csv", sep="."),  na="", col.names=TRUE, row.names=F, sep=",", quote=FALSE)

# load in Rqtl
library(qtl)
qtldata <- read.cross("csv", ".", "AfterJoinMapRemoving99SimCleaned20150331_0105.csv", genotypes=c("A","B"), na.strings=c(NA, "-"), estimate.map=F)
class(qtldata)[1] <-"riself"
# fix the order of chromosomes
qtldata$geno = qtldata$geno[c(4:14,1:3)]
summary(qtldata)

# now do some filtering that we couldn't do before without a rough genetic map..
# use "AfterGenMapFiltering.R"