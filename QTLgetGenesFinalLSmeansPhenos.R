# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

load("crossFinalLSmeansOrderFixedReRunCalcGeno_20151120_HIyield.Robject")
summary(crossH)
cross = crossH
load("modeldfAllFinalPhenosWithGrainHI_20151127_1630.RObject")
lsmeansmodeldf = modeldf
# remove epistatic QTL
lsmeansmodeldf = lsmeansmodeldf[!is.na(lsmeansmodeldf$chromosome),]
# change from factors to numbers
lsmeansmodeldf$chromosome = as.character(lsmeansmodeldf$chromosome)
lsmeansmodeldf$lowCIpos = as.numeric(as.character(lsmeansmodeldf$lowCIpos))
lsmeansmodeldf$hiCIpos = as.numeric(as.character(lsmeansmodeldf$hiCIpos))
# hmm maybe a way to loop through phenotypes and reps.. otherwise get an error about dup rows
# find.flanking(cross, chr=modeldfHTP$chromosome[1], pos=modeldfHTP$lowCIpos[1])
# modeldfHTP = modeldfHTP[modeldfHTP$phenotype == "NDVI" & modeldfHTP$rep == 1,]

# get the flanking for each CI
lowCIflanking = find.flanking(cross, chr=lsmeansmodeldf$chromosome, pos=lsmeansmodeldf$lowCIpos)
hiCIflanking = find.flanking(cross, chr=lsmeansmodeldf$chromosome, pos=lsmeansmodeldf$hiCIpos)

# add it to the table.. use right for low and left for hi
lsmeansmodeldf$lowCIflanking = lowCIflanking[,1]
lsmeansmodeldf$hiCIflanking = hiCIflanking[,2]

# if the CImarker is a real marker, keep it (because the flanking might be the next real marker)
# start by assumming we don't need to do anything
lsmeansmodeldf$lowCIflankingfinal = lsmeansmodeldf$lowCIflanking
lsmeansmodeldf$hiCIflankingfinal = lsmeansmodeldf$hiCIflanking

# make characters
as.data.frame(colnames(lsmeansmodeldf))
lsmeansmodeldf[ , c(14,15,20,21)] <- sapply(lsmeansmodeldf[ , c(14,15,20,21)], as.character)
            

library(gdata)

lsmeansmodeldf$lowCIflankingfinal[startsWith(lsmeansmodeldf$lowCImarker, "S")] = 
  lsmeansmodeldf$lowCImarker[startsWith(lsmeansmodeldf$lowCImarker, "S")]

lsmeansmodeldf$hiCIflankingfinal[startsWith(lsmeansmodeldf$hiCImarker, "S")] = 
  lsmeansmodeldf$hiCImarker[startsWith(lsmeansmodeldf$hiCImarker, "S")]

lsmeansmodeldf$lowCIflankingBP = 
  as.numeric(sapply(strsplit( lsmeansmodeldf$lowCIflankingfinal, "_"), "[[", 2))

lsmeansmodeldf$hiCIflankingBP = 
  as.numeric(sapply(strsplit( lsmeansmodeldf$hiCIflankingfinal, "_"), "[[", 2))

# give QTL names
# format: q then 2-5 cap letter for trait then chr ID then period then unique ID (so use peak pos)
lsmeansmodeldf$PhenoID = sapply(strsplit( as.character(lsmeansmodeldf$phenotype) , "_"), "[[", 1)
lsmeansmodeldf$PhenoID[lsmeansmodeldf$PhenoID == "biomass"] = "BMS"
lsmeansmodeldf$PhenoID[lsmeansmodeldf$PhenoID == "height"] = "PHT"
# lsmeansmodeldf$PhenoID[lsmeansmodeldf$PhenoID == "DTF"] = "FDN" # not sure about proper abbrev

lsmeansmodeldf$QTLname = paste0("q", lsmeansmodeldf$PhenoID, lsmeansmodeldf$chromosome, ".", round(lsmeansmodeldf$position))

# save it again
lsmeansmodeldfFINAL = lsmeansmodeldf[, c(26,2:4,7,8,12,23,24)]
filename = addStampToFilename("QTLsummaryForHei", "csv")
filename = addStampToFilename("FlankingBPForChr1QTLformanuscript", "csv")

write.table(lsmeansmodeldf, file=filename, col.names=T, row.names=F, sep=",", quote=FALSE)

# ### get genes #####
# ###################
# 
# # first get gff file with genes from MSU..
# fileurl = "ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.gff3"
# download.file(fileurl, "msu7.gff3")
# setwd("~/Desktop/Dropbox/working folder/Rprojects/RILpop/data/")
msu7gff = read.delim("msu7.gff3", header=F, comment.char="#")
# # proper names
colnames(msu7gff) = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
msu7gff$seqid = as.character(msu7gff$seqid)
# # maybe a faster way to do this?
split = strsplit(msu7gff$seqid, "Chr")
msu7gff$chr = sapply(split, "[", 2)
# # only keep genes
msu7gff = msu7gff[msu7gff$type == "gene",]
# # and get the gene ID
split = strsplit(as.character(msu7gff$attributes), ";")
split2 = sapply(split, "[", 1)
split3 = strsplit(split2, "ID=")
msu7gff$gene.id = sapply(split3, "[", 2)

# get genes
# as.data.frame(colnames(modeldfHTP))
# i = 1
# gff.in <- msu7gff[msu7gff$chr == modeldfHTP[i,"chromosome"] & 
#                     msu7gff$end <= modeldfHTP[i,"hiCIflankingBP"] & 
#                     msu7gff$start >= modeldfHTP[i,"lowCIflankingBP"] , ]

# modeldfHTP[1,c(5,6,24,25)]
# summary(gff.in)

# this takes a while.. maybe faster with a matrix or something
genelist = data.frame()
# for each QTL

# just for Chr 3
lsmeansmodeldf_chr3 = lsmeansmodeldf[37,]
lsmeansmodeldfbackup = lsmeansmodeldf
lsmeansmodeldf = lsmeansmodeldf_chr3

for (i in 1: length(rownames(lsmeansmodeldf))){
  # put list of genes in df
  # where the chr are the same, and range of bp match
  gff.in <- msu7gff[msu7gff$chr == lsmeansmodeldf[i,"chromosome"] & 
                      msu7gff$end <= lsmeansmodeldf[i,"hiCIflankingBP"] & 
                      msu7gff$start >= lsmeansmodeldf[i,"lowCIflankingBP"] , ]
  # if empty, skip
  if(dim(gff.in)[1] == 0){print ("skipping");next}
  # otherwise add QTL data
  gff.in = cbind(gff.in, lsmeansmodeldf[i, c(24, 2:3,6,7,11, 22, 23)])
  # and combine into big df
  genelist = rbind(genelist, gff.in)
}

# save it
filename = addStampToFilename("CompleteGeneListChr3_QTL_ForHei", "tsv")
write.table(genelist, file=filename, col.names=T, row.names=F, sep="\t", quote=FALSE)
filename = addStampToFilename("CompleteGeneListFinalLSmeansQTL", "Robject")
save(genelist, file=filename)

# filter and reorganize it
genelistfinal = genelist
as.data.frame(colnames(genelistfinal))
# genelistfinal = genelistfinal[,c(13,12,14:18,11,9)]
genelistfinal = genelistfinal[,c(12:19,11,9)]


# parse the attributes
genelistfinal$GO = sapply(strsplit(as.character(genelistfinal$attributes), ";Note="), "[[", 2)

# fix url encoding http://www.w3schools.com/tags/ref_urlencode.asp
genelistfinal$GOfixed = gsub("%20", " ", genelistfinal$GO )
genelistfinal$GOfixed = gsub("%2C", ",", genelistfinal$GOfixed )
genelistfinal$GOfixed = gsub("%2F", "/", genelistfinal$GOfixed )
genelistfinal$GOfixed = gsub("%2F", "*", genelistfinal$GOfixed )
genelistfinal$GOfixed = gsub("%2B", "+", genelistfinal$GOfixed )
genelistfinal$GOfixed = gsub("%3B", ";", genelistfinal$GOfixed )
genelistfinal$GOfixed = gsub("%3A", ":", genelistfinal$GOfixed )
genelistfinal$GOfixed = gsub("%2A", "*", genelistfinal$GOfixed )
genelistfinal$GOfinal = gsub("%27", "'", genelistfinal$GOfixed )

# genelistfinal = genelistfinal[,-c(9:11)]
genelistfinal = genelistfinal[,-c(10:12)]


filename = addStampToFilename("CompleteGeneList_Chr3_QTLFormatted", "tsv")
write.table(genelistfinal, file=filename, col.names=T, row.names=F, sep="\t", quote=FALSE)

filename = addStampToFilename("CompleteGeneListFinalLSmeansQTLFormatted", "Robject")
save(genelistfinal, file=filename)

# 
# # lovell's approach
# # read in the gff file to get genes within each region
# marker.info<-read.csv("marker.info.csv")
# mark.info<-marker.info
# 
# # mark.info = markerlist
# mark.info$qtl<-"no"
# statsbp.out<-data.frame()
# for(i in unique(beststats$qtl.id)){
#   stats<-beststats[beststats$qtl.id==i,]
#   mark.info.low<-mark.info[mark.info$chr==as.character(stats$chromosome) &
#                              mark.info$cm<=as.numeric(stats$lowCIpos),]
#   mark.info.low<-mark.info.low[which(abs(mark.info.low$cm-as.numeric(stats$lowCIpos)) == 
#                                        min(abs(mark.info.low$cm-as.numeric(stats$lowCIpos)))),]
#   mark.info.hi<-mark.info[mark.info$chr==as.character(stats$chromosome) &
#                             mark.info$cm>=as.numeric(stats$hiCIpos),]
#   mark.info.hi<-mark.info.hi[which(abs(mark.info.hi$cm-as.numeric(stats$hiCIpos)) == 
#                                      min(abs(mark.info.hi$cm-as.numeric(stats$hiCIpos)))),]
#   if(length(mark.info.hi[,1])>1){mark.info.hi<-mark.info.hi[mark.info.hi$bp==max(mark.info.hi$bp),]}
#   if(length(mark.info.low[,1])>1){mark.info.low<-mark.info.low[mark.info.low$bp==min(mark.info.low$bp),]}
#   #select all markers based on cm
#   marker.info.ci<-mark.info[mark.info$cm>=mark.info.low$cm & mark.info$cm<=mark.info.hi$cm & mark.info$chr==as.character(stats$chromosome),]
#   mark.info$qtl[mark.info$marker %in% marker.info.ci$marker]<-"QTL"
#   max.bp<-max(marker.info.ci$bp)
#   min.bp<-min(marker.info.ci$bp)
#   stats$lowbp<-min.bp; stats$hibp<-max.bp
#   statsbp.out<-rbind(statsbp.out, stats)
# }
# 
# # fix linkage group names to match
# statsbp.out$chr = as.numeric(substr(statsbp.out$chromosome, 1, 2))
# 
# # ggplot(mark.info, aes(x=bp, y=cm, color=qtl, fill=qtl,alpha=qtl))+
# #   geom_point()+
# #   scale_fill_manual(values=c("grey","red"))+
# #   scale_alpha_manual(values=c(1, 0.1,.1))+
# #   facet_wrap(~chr)
# 
# gff = msu7gff
# 
# for (i in 1: length(rownames(statsbp.out))){
#   gff.in<-gff[gff$chr==statsbp.out[i,"chr"] & gff$end<=statsbp.out[i,"hibp"] & gff$start>=statsbp.out[i,"lowbp"],]
#   #   gff.in$qtl = statsbp.out[i,18]
#   genelist = rbind(genelist, gff.in)
#   candidate.genes<-as.character(gff.in$gene.id)
#   #   print(candidate.genes)
# }
# 
# # not sure why duplicates
# finalgenelist = unique(genelist)
# 
# fileprefix = paste("genelist_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(finalgenelist, file=paste(fileprefix, "tsv", sep="."),row.names=F, sep="\t", quote=F)
