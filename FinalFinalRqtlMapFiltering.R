# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

load("qtldataAfterRemoveBadChr12NewJoinMapPreFilter_20150506_1701.Robject")

qtldataB4filter = qtldata

load("qtldataAfterRemoveBadChr12NewJoinMapB4ripand1markerpercM_20150506_1713.Robject")
qtldataAfterSegDist = qtldata

pdf("geneticmaps_B4andAfterSegDistFiltering.pdf", width=10, height=8)
plot.map(qtldataB4filter, qtldataAfterSegDist, main="B4 filtering on left, after seg dist filtering on right")
dev.off()

qtldata = qtldataB4filter

# plot seg dist
plot.map(qtldata)
map<-pull.map(qtldata, as.table=T)
map$gene<-rownames(map)
gt<-geno.table(qtldata)
gt$gene<-rownames(gt)
gt$prop.AA<-gt$AA/(gt$AA+gt$BB)
gt$prop.BB<-gt$BB/(gt$AA+gt$BB)
gt$gene<-rownames(gt)
gt$log10p<-log10(gt$P.value)
gt<-merge(gt,map, by=c("gene","chr"))

segdistplot = ggplot(gt, aes(x=pos,y=log10p)) +
  geom_point() + geom_line() + facet_wrap(~chr) +
  scale_y_reverse("log10 chi2 p.value for seg. dist.") +
  scale_x_continuous("mapping position (cM)") +
  theme_bw() + 
  ggtitle("distribution of segregation distortion in the map log10 p values")

segdistplot

# plot ratios instead..
segdistratiosplot = ggplot(gt, aes(x=pos)) +
  geom_line(aes(y = prop.AA), color="red") +
  geom_point(aes(y = prop.AA), color="red") + 
  geom_line(aes(y = prop.BB), color="blue") +
  geom_point(aes(y = prop.BB), color="blue") + 
  facet_wrap(~chr) +
  scale_y_continuous("%") +
  scale_x_continuous("mapping position (cM)") +
  theme_bw() + 
  ggtitle("distribution of segregation distortion in the map Aswina=red IR64=blue")

segdistratiosplot

ggsave("segdistratios_B4filtering.pdf", segdistratiosplot)


# subset markers, keep ones based on least seg dist and least missing data..
# I think this is my final list of markers?

out<-vector()
for(i in chrnames(qtldata)){
  gt<-geno.table(qtldata, chr=i)
  n.missing <- nmissing(subset(qtldata, chr=i), what="mar")
#   sd<- scale(gt$P.value) # don't do this.. it removes too many markers...
#   miss <- scale(-log( (n.missing+1) / (nind(qtldata)+1) ))
  wts <- -log( (n.missing+1) / (nind(qtldata)+1) )
#   wts<-sd+miss
  ms<-pickMarkerSubset(pull.map(qtldata)[[i]], .5, wts)
  out<-c(out, ms)
}

# so I need to take this list and make this my final list of markers?
markers = row.names(map)
tokeep = markers[markers %in% out]
todiscard = markers[!(markers %in% out)]
qtldataAfterPickMarkerSubsetNoSD = drop.markers(qtldata, todiscard)
summary(qtldataAfterPickMarkerSubsetNoSD)

# compare..
pdf("geneticmaps_B4andAfterSegDistFilteringWithoutSD.pdf", width=10, height=8)
plot.map(qtldataB4filter, qtldataAfterPickMarkerSubsetNoSD, main="B4 filtering on left, after seg dist w/o SD on right")
dev.off()

# now redo the seg dist plot...
plot.map(qtldataAfterPickMarkerSubsetNoSD)
map<-pull.map(qtldataAfterPickMarkerSubset, as.table=T)
map$gene<-rownames(map)
gt<-geno.table(qtldataAfterPickMarkerSubset)
gt$gene<-rownames(gt)
gt$prop.AA<-gt$AA/(gt$AA+gt$BB)
gt$prop.BB<-gt$BB/(gt$AA+gt$BB)
gt$gene<-rownames(gt)
gt$log10p<-log10(gt$P.value)
gt<-merge(gt,map, by=c("gene","chr"))

segdistratiosplot = ggplot(gt, aes(x=pos)) +
  geom_line(aes(y = prop.AA), color="red") +
  geom_point(aes(y = prop.AA), color="red") + 
  geom_line(aes(y = prop.BB), color="blue") +
  geom_point(aes(y = prop.BB), color="blue") + 
  facet_wrap(~chr) +
  scale_y_continuous("%") +
  scale_x_continuous("mapping position (cM)") +
  theme_bw() + 
  ggtitle("seg dist after pickmarkersubset Aswina=red IR64=blue")

segdistratiosplot
ggsave("segdistratios_afterfiltering.pdf", segdistratiosplot)

# check markers at end of chr3..
subset(map, chr=="3")

# overwrite qtldata now that I think this worked
qtldata = qtldataAfterPickMarkerSubsetNoSD

filename = addStampToFilename("qtldataAfterRemoveBadChr12NewJoinMapB4ripand1markerpercM", "Robject")
save(qtldata, file=filename)

load("qtldataAfterRemoveBadChr12NewJoinMapB4ripand1markerpercM_20150506_1713.Robject")

# check marker order within a window
# start with countxo over large window, 
# then repeat with smaller window with likelihood (slower)
cross.reorder<-qtldata
rips<-list()
# this takes a minute or so..
for(i in chrnames(cross.reorder)){
  rips[[i]]<-ripple(cross.reorder, chr=i, window=5, method="countxo", error.prob=0.001,
                    map.function="kosambi", verbose=T)
}

# for each chr, pull the original, and best order
ripchr = data.frame()
for(i in chrnames(cross.reorder)){
  temp = as.data.frame(rips[[i]][c(1:2), ncol(rips[[i]])])
  temp$chr = i
  ripchr = rbind(ripchr, temp)
}
# fix columns
ripchr = as.data.frame(cbind(order = rownames(ripchr), ripchr)) 
names(ripchr) = c("order", "obligoXO", "chr")

# export
# filename = addStampToFilename("RippleResults", "tsv")
# setwd(datapath)
# write.table(ripchr, file=filename , col.names=TRUE, row.names=F, sep="\t", quote=FALSE)

# figure out when original is not best order
# make it wide to check min and add to list of chr
library(reshape2)
ripchrwide = dcast(ripchr, order ~ chr, value.var="obligoXO")
head(ripchrwide)
# need this for starts with
library(gdata)
# skip order col
i = 2
toswitch = vector()
# for each col
for (i in i:length(colnames(ripchrwide))){
  # find row index for min value for each col
  rowidx = which(ripchrwide[,i] == min(ripchrwide[,i], na.rm=T))
  # add to list to switch order unless Initial is min, or there is a tie
  if(!startsWith(ripchrwide[rowidx,1], "Initial",) && length(rowidx) < 2){
    toswitch = c(toswitch, colnames(ripchrwide)[i])
  }
}

# make the switches
i=1
for(i in toswitch){
  cross.reorder <- switch.order(cross.reorder, i, rips[[i]][2,])
}


qtldata = cross.reorder
setwd(datapath)
filename = addStampToFilename("qtldataAfterRipB41cM", "Robject")
save(qtldata, file=filename)
load("qtldataAfterRipB41cM_20150507_1741.Robject")
qtldata = est.rf(qtldata)
filename = addStampToFilename("EstRFmapbefore1cMmarkerAfterRemoveBadChr12", "png")
png(filename, width=9, height=8.5, units="in", res=150)
# warning about chr is ok..
plotRF(qtldata, par(mar=c(5.1, 4.1, 4.1, 4)))
dev.off()

# check map.. not sure I should do this now..
qtldata2 = est.rf(qtldata)
# pdf suck for these maps
# make margins big enough for names
par()$mar
# bottom, left, top, and right
# par(mar=c(5.1, 4.1, 4.1, 4))
filename = addStampToFilename("EstRFmap", "png")
png(filename, width=9, height=8.5, units="in", res=150)
# warning about chr is ok..
plotRF(qtldata2, par(mar=c(5.1, 4.1, 4.1, 4)))
dev.off()


# find the markers from the plot of the estRF...
temp<-pull.rf(qtldata2, what="lod", chr=12)
# the entire chr
plot(temp[1:72])
# the bad part
plot(temp[1:25])

# then pull these bad markers from the map.. and redo the map in joinmap..

# pull entire rf matrix, pull sub matrix between chrs of interest, plot heatmap, and identify regions to remove
chr11 = markernames(qtldata2, chr=11)
chr12 = markernames(qtldata2, chr=12)
temp2 = pull.rf(qtldata2)
asdf = temp2[chr11,chr12]
par(mar=c(5.1, 4.1, 4.1, 4))
par( "mar" ) + c( 2, 4, 0, 0 ) )
image(asdf, xaxt = "n", yaxt = "n")
# las is label orientation
axis( 1, at=seq(0,1,length.out=ncol( asdf ) ), labels= colnames( asdf ), las= 2 , cex=0.7)
axis( 2, at=seq(0,1,length.out=nrow( asdf ) ), labels= rownames( asdf ), las= 2 , cex=0.8)

# id linked markers

# first get col index (I suppose you could graph it with index above..)
colidx = grep("^S12_23919524$", colnames(asdf))
linkedmarkers = colnames(asdf)[colidx:(colidx+12)]

image(asdf, xaxt = "n", yaxt = "n")
# las is label orientation
axis( 1, at=seq(0,1,length.out=ncol( asdf ) ), labels= col( asdf )[1,], las= 2 , cex=0.7)
axis( 2, at=seq(0,1,length.out=nrow( asdf ) ), labels= rev(row( asdf )[,1]), las= 2 , cex=0.8)

linkedmarkers = rownames(asdf)[1:12]
linkedmarkers = c(linkedmarkers, colnames(asdf)[31:45])

# check these
qtldataDropLinked = drop.markers(qtldata2, linkedmarkers)
summary(qtldataDropLinked)
summary(qtldata2)
qtldata2 = qtldataDropLinked
#####################
# not sure this works or is necessary..
#####################

cross = qtldata
# keep one marker per cM based on which markers have least missing data
out<-vector() 
for(i in chrnames(cross)){
  gt<-geno.table(cross, chr=i)
  n.missing <- nmissing(subset(cross, chr=i), what="mar")
  wts <- -log( (n.missing+1) / (nind(cross)+1) )
  ms<-pickMarkerSubset(pull.map(cross)[[i]], 1, wts) # one marker per cM..
  out<-c(out, ms)
}
cull.cross<-drop.markers(cross, markernames(cross)[!markernames(cross) %in% out])
# newmap<-est.map(cull.cross, error.prob=0.001, map.function="kosambi")
filename = addStampToFilename("mapsB4andAfterEst.Map", "pdf")
pdf(filename, width=10, height=8)
plot.map(cull.cross, newmap, main="original map left, after est.map right")
dev.off()
# not sure I want to use the map that R generates.. it is smaller.. no other difference..
# cull.cross <- replace.map(cull.cross, newmap)
filename = addStampToFilename("qtldataFinalWithRestmap", "Robject")
save(cull.cross, file=filename)
filename = addStampToFilename("mapsB4andAfter1markerpercM", "pdf")
setwd(datapath)
pdf(filename, width=10, height=8)
plot.map(cross, cull.cross, main="original map left, after 1 marker per cM right")
dev.off()
cross<-cull.cross

setwd(datapath)
filename = addStampToFilename("qtldataAfterSegDistnoSDRip1markerpercM", "Robject")
save(cross, file=filename)
load("qtldataAfterSegDistRip1markerpercM_20150429_1525.Robject")
cross = est.rf(cross)
filename = addStampToFilename("EstRFmapFinal", "png")
png(filename, width=9, height=8.5, units="in", res=150)
# warning about chr is ok..
plotRF(cross, par(mar=c(5.1, 4.1, 4.1, 4)))
dev.off()