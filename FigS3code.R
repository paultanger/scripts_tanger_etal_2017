# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

rm(list=ls())
load("modeldfHTP_20150616_1549.Robject")
load("crossFinalLSmeansOrderFixedReRunCalcGeno_20161002_newHI.Robject")

library(qtl)
phes<-c("DTF_LSmeans", "DTH_LSmeans", "biomass_LSmeans", "height_LSmeans", "grain_rep3", "HI_rep3")
s1s<-scanone(crossH, method="hk", phe=phes)

# top of Figure S3
pdf("newHIScanonePlots_allchr_regyscale.pdf", height=6, width=6)
plot(s1s, bty="n", ylim=c(0,160), ylab="LOD Score", type="n")
palette<-c("black","pink","darkorange","cyan","darkgreen","blue")
library(RColorBrewer)
pal<-brewer.pal(name="Paired", n=6)
for(i in 1:6) plot(s1s, lodcolumn=i, add=T, col=pal[i])
rect(xleft=0, xright=1900, ybottom=0, ytop=3.32, col=rgb(0,0,0,.5))
#abline(h=log10(3.32+1), col="grey", lty=2)
legend("topright",phes, col=pal, lty=1)
dev.off()

htp<-modeldfHTP[,c("phenotype", "rep", "chromosome", "position")]
htp2<-htp[htp$phenotype %in% c("Chla","CTD","HTPheight","NDRE","NDVI"),]

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

# bottom of Figure S3
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
