# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")
rm(list=ls())
load("modeldfHTP_20150616_1549.Robject")
load("crossFinalLSmeansOrderFixedReRunCalcGeno_20151120_HIyield.Robject")

phes<-c("DTF_LSmeans", "DTH_LSmeans", "biomass_LSmeans", "height_LSmeans", "grain_rep3", "HI_rep3")
s1s<-scanone(crossH, method="hk", phe=phes)

pdf("newScanonePlots_chr136_regyscale.pdf", height=3, width=6)
plot(s1s, bty="n", ylim=c(0,160), ylab="LOD Score", type="n", chr=c(1,3,6))
palette<-c("black","pink","darkorange","cyan","darkgreen","blue")
library(RColorBrewer)
pal<-brewer.pal(name="Paired", n=6)
for(i in 1:6) plot(s1s, lodcolumn=i, add=T, col=pal[i], chr=c(1,3,6))
abline(h=log10(3.32+1), col="grey", lty=2)
#legend("topright",phes, col=pal, lty=1)
dev.off()

pdf("newScanonePlots_allchr_regyscale.pdf", height=6, width=6)
plot(s1s, bty="n", ylim=c(0,160), ylab="LOD Score", type="n")
palette<-c("black","pink","darkorange","cyan","darkgreen","blue")
library(RColorBrewer)
pal<-brewer.pal(name="Paired", n=6)
for(i in 1:6) plot(s1s, lodcolumn=i, add=T, col=pal[i])
rect(xleft=0, xright=1900, ybottom=0, ytop=3.32, col=rgb(0,0,0,.5))
#abline(h=log10(3.32+1), col="grey", lty=2)

#legend("topright",phes, col=pal, lty=1)
dev.off()



s1s.log10<-s1s
for(i in phes) s1s.log10[,i]<-log10(s1s.log10[,i]+1)
maxy<-log10(250+1)
pdf("newScanonePlots_chr136.pdf", height=3, width=6)
plot(s1s.log10, yaxt="n", bty="n", ylim=c(0,maxy), ylab="LOD Score (log10 scale)", type="n", chr=c(1,3,6))
axis(2, at=c(0, log10(10+1), log10(50+1), log10(100+1),maxy), labels=c(0,10,50,100,200))
palette<-c("black","pink","darkorange","cyan","darkgreen","blue")
library(RColorBrewer)
pal<-brewer.pal(name="Paired", n=6)
for(i in 1:6) plot(s1s.log10, lodcolumn=i, add=T, col=pal[i], chr=c(1,3,6))
abline(h=log10(3.32+1), col="grey", lty=2)
legend("topleft",phes, col=pal, lty=1)
dev.off()

htp<-modeldfHTP[,c("phenotype", "rep", "chromosome", "position")]
htp<-htp[htp$chromosome %in% c(1,3,6),]
htp$chromosome<-factor(htp$chromosome, levels=c(1,3,6))

htp2<-htp[htp$phenotype %in% c("Chla","CTD","HTPheight","NDRE","NDVI"),]
phes.htp<-unique(htp2$phenotype)
out.list<-list()
tos1<-lapply(phes.htp, function(j){
  outall<-data.frame()
  for(i in c(1,3,6)){
    dat<-htp2[htp2$chromosome==i & htp2$phenotype==j,]
    if(nrow(dat)>0){
      end<-max(dat$position)
      out<-hist(dat$position, breaks=seq(from=0, to=ceiling(end/4)*4, by=4))
      out<-data.frame(chr=i, pos=out$mids, phe=out$counts)
      colnames(out)[3]<-paste("count",j,sep="_")
      outall<-rbind(outall,out)
    }
  }
  return(outall)
})
library(plyr)
library(reshape)
tos2<-merge_recurse(tos1, by=c("chr","pos"))
test<-merge(s1s.log10, tos2)

class(tos2)<-c("scanone","data.frame")
map<-pull.map(crossH, as.table=T)
colnames(map)<-c("chromosome","position")
map<-map[map$chromosome %in% c(1,3,6),]
map$map=map$position
map<-map[,-2]
map$rep<-NA
map$phenotype<-NA
map$position<-NA
htp2$map<-NA
map<-map[,colnames(htp2)]
htp3<-rbind(htp2,map)

htp3<-merge(htp2, map, by=c("chromosome"),all.x=T)
pdf("HTP_QTL_histograms_coloredByRep.pdf",height=4, width=6)
ggplot(htp3, aes(x=position, fill=rep))+
  scale_fill_manual(values=c("pink","cornflowerblue","darkred"))+
  geom_histogram(binwidth=4,origin = 0)+
  facet_grid(phenotype~chromosome, scales="free", space="free_x")+
  theme_bw()+
  geom_rug(sides="b", aes(x=map))+
  theme(strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
dev.off()

pdf("HTP_QTL_histograms_coloredByPhe.pdf",height=2, width=6)
ggplot(htp3, aes(x=position, fill=phenotype))+
  scale_fill_manual(values=pal)+
  geom_histogram(binwidth=4,origin = 0)+
  facet_grid(~chromosome, scales="free", space="free_x")+
  theme_bw()+
  geom_rug(sides="b", aes(x=map))+
  theme(strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
dev.off()


htp<-modeldfHTP[,c("phenotype", "rep", "chromosome", "position")]
htp2<-htp[htp$phenotype %in% c("Chla","CTD","HTPheight","NDRE","NDVI"),]

class(tos2)<-c("scanone","data.frame")
map<-pull.map(crossH, as.table=T)
colnames(map)<-c("chromosome","position")
map$map=map$position
map<-map[,-2]
map$rep<-NA
map$phenotype<-NA
map$position<-NA
htp2$map<-NA
map<-map[,colnames(htp2)]
htp3<-rbind(htp2,map)
htp3<-htp3[!is.na(htp3$chromosome),]
pdf("HTP_QTL_histograms_coloredByRep_allchr.pdf",height=4, width=6)
htp3$chromosome<-as.numeric(as.character(htp3$chromosome))
ggplot(htp3, aes(x=position, fill=rep))+
  scale_fill_manual(values=c("pink","cornflowerblue","darkred"))+
  geom_histogram(binwidth=10,origin = 0)+
  facet_grid(phenotype~chromosome, scales="free", space="free_x")+
  theme_bw()+
  geom_rug(sides="b", aes(x=map))+
  theme(strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
dev.off()