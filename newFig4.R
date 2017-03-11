# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

tp<-tp[complete.cases(tp),]
tp2<-tp[tp$phenotype %in% c("HI_rep3","grain_rep3") & complete.cases(tp),]
tp2$unique<-paste(tp2$chromosome, tp2$position, sep="_")

tp.grain<-tp[tp$phenotype == "grain_rep3",]
tp.hi<-tp[tp$phenotype == "HI_rep3",]
names(tp.hi)<-paste(names(tp.hi),"hi", sep="_")
names(tp.grain)<-paste(names(tp.grain),"grain", sep="_")
library(ggplot2)
colnames(tp.hi)[2:3]<-c("chr","pos")
colnames(tp.grain)[2:3]<-c("chr","pos")
tp2<-merge(tp.hi,tp.grain, by=c("chr","pos"),all.x=T)
ggplot(tp2, aes(x=chromosome, y=position))+geom_point()


tp.grain<-tp[tp$phenotype == "grain_rep3",]
tp.hi<-tp[tp$phenotype == "HI_rep3",]
m.hi<-makeqtl(cross, chr=tp.hi$chromosome, pos=tp.hi$position, what="prob")
m.grain<-makeqtl(cross, chr=tp.grain$chromosome, pos=tp.grain$position, what="prob")
ref.hi<-refineqtl(cross, qtl=m.hi, pheno.col="grain_3_final", method="hk")
plotLodProfile(ref.hi)

s.grain<-stepwiseqtl(cross, pheno.col="grain_3_final", method="hk", additive.only=T, penalties=calc.penalties(perms.out[[6]]), max.qtl=8)
s.hi<-stepwiseqtl(cross, pheno.col="HI_3_final", method="hk", additive.only=T, penalties=calc.penalties(perms.out[[7]]), max.qtl=8)
par(mfrow=c(2,1))
plotLodProfile(s.grain)
plotLodProfile(s.hi)
cons.mod<-makeqtl(cross, chr=c(1,1,2,3,4,4,5,6,6,9),
                  pos=c(26,166,51,10,89,109,35,22.5,52,19), what="prob")
fit.hi<-fitqtl(cross, qtl=cons.mod, method="hk", pheno.col="HI_3_final", get.ests=T)
fit.grain<-fitqtl(cross, qtl=cons.mod, method="hk", pheno.col="grain_3_final", get.ests=T)
ref.grain<-refineqtl(cross, qtl=m.hi, pheno.col="grain_3_final", method="hk")
est.grain<-data.frame(summary(fit.grain)$ests)
est.hi<-data.frame(summary(fit.hi)$ests)
est.grain$qtl<-row.names(est.grain)
est.hi$qtl<-row.names(est.hi)
colnames(est.hi)<-paste("hi",colnames(est.hi), sep="_")
colnames(est.grain)<-paste("grain",colnames(est.grain), sep="_")
ests<-cbind(est.grain, est.hi)
ests<-ests[!rownames(ests)=="Intercept",]
plot(ests$hi_est, ests$grain_est, bty="n")
abline(h=0)
abline(v=0)
est.hi$phe<-"hi"
est.grain$phe<-"grain"

ndvi.phes<-phenames(cross)[grep("NDVI_1", phenames(cross))]
fits<-lapply(ndvi.phes, function(x){
  fit<-fitqtl(cross, qtl=cons.mod, method="hk", pheno.col=x, get.ests=T)
  est<-data.frame(summary(fit)$ests)
  est$qtl<-row.names(est)
  est<-est[!rownames(est)=="Intercept",]
  return(est)
})
names(fits)<-ndvi.phes
fits2<-ldply(fits, data.frame)
ests$ndvi.eff<-tapply(fits2$est, fits2$qtl, mean)
ggplot(ests, aes(x=hi_est, y=grain_est, col=ndvi.eff))+geom_point(size=5)+
  theme_classic()+
  geom_vline(xintercept = 0)+
  geom_hline(xintercept = 0)+
  theme(legend.justification = c(0, 1), legend.position = c(0, 1),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  scale_colour_gradient(low="black",high="red")
data.frame(ests)
rbind(est.hi, est.grain)
plotLodProfile(ref.grain)

plot(s1s.2, lodcolumn=7:9, type="h")

hi.mod<-makeqtl(cross, chr=tp.hi)
for(j in poss){
  temp<-dfdat[dfdat$lowCIpos<=j & dfdat$hiCIpos>=j,]
  
}
}

for(i in 1:12){
  tp$position[tp$chr==i]<-tp$position[tp$chr==i]+c(0,cumsum(chrlen(cross)+25)[1:8])[i]
}

nbins<-nrow(s1s)
counts<-hist(tp$position, breaks=nbins, plot=F)