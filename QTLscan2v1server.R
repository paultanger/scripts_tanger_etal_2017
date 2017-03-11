# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

# this is for scan2...

load("allphenosscan2perms_20150603_2334.Robject")
load("phenoswithscan1QTLalpha0.1_20150528_1728.Robject")
load("crossallphenosOrderFixedReRunCalcGeno_20150612_1559.Robject")
cross = cross2
# lovell
setwd(sourcepath)
setwd("../r-QTL_functions/")
setwd("../lovells_functions_githubfork/")
source("SourceAll.R")
sourceDir("../r-QTL_functions/")
sourceDir("../lovells_functions_githubfork/")
setwd(datapath)

# so perms.out is a list with perms for each phenotype..
# pens = calc.penalties(perms.out[[1]], alpha=0.05)
# pens

# do it for all
# pens = lapply(perms.out, function(x) calc.penalties(x, alpha=0.05))
# pens[1]

# lovell's method
pens.out<-data.frame()
for (i in qtlphes){pens<-calc.penalties(perms.out[[i]], alpha=.1); pens.out<-rbind(pens.out,pens)}
rownames(pens.out)<-qtlphes; colnames(pens.out)<-c("main","heavy","light")
pens.out

maxqtl = 10
# maxqtl = 6

modelList<-list()
#run model selection

# for (i in qtlphes){
#   stepout<-stepwiseqtl(cross, pheno.col=i, method="hk",
# #                                       penalties=pens[i], 
#                                       penalties= pens.out[i,],
#                                       max.qtl=maxqtl, 
#                                       covar=NULL,
#                                       refine.locations=T,
#                                       additive.only=F, 
#                                       keeplodprofile=T,
#                                       keeptrace=T,
#                                       verbose=F)
#   models.out[[i]]<-stepout
# }
nodes = 23
# or do it parallelized.. errors out on server with 3.0.2 but works for 3.1.2
ptm <- proc.time()
modelList <- mclapply(qtlphes, mc.cores=nodes, function (i) { # could try 
  stepout <- stepwiseqtl(cross, pheno.col=i, method="hk",
                       penalties=pens.out[i,],
                       max.qtl=maxqtl,
                       covar=NULL,
                       refine.locations=T,additive.only=F, 
                       keeplodprofile=T,keeptrace=T,verbose=F)
})
proc.time() - ptm

# need to name elements in list
names(modelList) = qtlphes

filename = addStampToFilename("stepwiseoutputMaxQTL10AfterFixChrOrder", "Robject")
save(modelList, file=filename)

# btw, to check perms for a specific pheno:
summary(perms.out[["NDVI_3_030813"]], alpha=0.05)

load("stepwiseoutputMaxQTL10AfterFixChrOrder_20150612_1722.Robject")


# get data into a df
stepParsed<-lapply(qtlphes, function (i) {
  stepwiseqtl.summary.simple(cross, model.in= modelList[[i]], phe=i, covar=NULL, ci.method="drop", drop=1.5, plot=FALSE, printout=T)
})
stepParsed

# try new function
newmodel = lapply(qtlphes, function (i) {
  stepwiseStats(cross, model.in= modelList[[i]], phe=i, covar=NULL, ci.method="drop", drop=1.5, plot=FALSE, printout=T)
})
newmodeldf <- ldply(newmodel, data.frame)
# doesn't work right now
library(plyr)
modeldf <- ldply(stepParsed, data.frame)

setwd(datapath)
filename = addStampToFilename("stepwiseDFafterfixchr", "RObject")
save(modeldf, file=filename)
filename = addStampToFilename("stepwiseDFafterfixchr", "csv")
write.table(modeldf, file=filename, col.names=T, row.names=F, sep=",", quote=FALSE)

# do this with bayes CI
stepParsedwithBayes <- lapply(qtlphes, function (i) {
  stepwiseqtl.summary.simple(cross, model.in= modelList[[i]], phe=i, covar=NULL, ci.method="bayes", drop=1.5, plot=FALSE, printout=T, prob=.9)
})
stepParsedwithBayes
library(plyr)
modeldfwithBayes <- ldply(stepParsedwithBayes, data.frame)

filename = addStampToFilename("stepwiseDFwithBayes", "csv")
write.table(modeldfwithBayes, file=filename, col.names=T, row.names=F, sep=",", quote=FALSE)

load("stepwiseDFafterfixchr_20150616_1213.RObject")
# try lovell's plot function
library(ggplot2)
fullqtlplotnobayes = plot.qtlci(  cross, # cross object
             modeldf, # df of stats
             "phenotype", # these are names of cols in df for phenotypes
             "chromosome", # these are names of cols in df for chr
             "lowCIpos", # these are names of cols in df for left CI
             "position", # these are names of cols in df for center pos
             "hiCIpos", # these are names of cols in df for right CI
                       models=modelList, # list of models
                       plottype="lines", # lines or lodprofile
                       standardize=TRUE,
                       ci.thickness=1.5,
                       point.size=2
             )

filename = addStampToFilename("fullQTLlinesplot", "pdf")
ggsave(filename, fullqtlplot, width=25, height=20, units="in", scale=2, limitsize=F)


# new function
# plotMultiQTL(cross=cross, stats=modeldf, pointsize=.5, pointshape=FALSE, linethickness=2, colbychr=T, palette=rainbow, background=TRUE)
library(plyr)
library(gdata)
phenostoplot = modeldf[startsWith(modeldf$phenotype, pattern ="NDVI"),]
phenostoplot = subset(phenostoplot, chromosome == 2)

filename = addStampToFilename("chr2looks weird", "pdf")
pdf(filename, width=25, height=20)
newqtlplotnobayes = plotMultiQTL(  cross, # cross object
                                   phenostoplot, # df of stats
                                  pointsize=.5, 
                                  pointshape=1, 
                                  linethickness=2, 
                                  colbychr=T, 
                                  palette=rainbow, 
                                  background=TRUE
)
dev.off()

# try with new function
library(gdata)
library(plyr)
# phenostoplot = levels(modeldf$phenotype)[startsWith(levels(modeldf$phenotype), pattern ="NDVI")]
# subset df before..

# what if we use DAS instead?
load("modeldfHTP_20150616_1549.Robject")
phenostoplot = modeldfHTP[startsWith(modeldfHTP$phenotype, pattern ="NDVI"),]
# rename col to trick plot function to use DAS on y axis
phenostoplot = phenostoplot[,-2]
colnames(phenostoplot)[20] = "phenotype"
# sort by DAS..
phenostoplot = phenostoplot[with(phenostoplot, order(-phenotype)),]

# create colors for LOD score
# install.packages("RColorBrewer")
library(RColorBrewer)
# how many colors do you want?
numcolors = 7
redscale = colorRampPalette(brewer.pal(6,"Reds"))(numcolors)
# the lowest color is too faint
mypalette = redscale[2:numcolors]
# two colors
# mypalette = c("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b")
# one color
# mypalette = c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#b30000")

# divide LOD into rounded equal groups..
# add 1 to max since we are rounding (prevents NAs)
groups = round(seq(0, round(max(phenostoplot$LOD))+1, length.out = numcolors))

# create color vector based on num of colors you want and LOD scores
palettevector = mypalette[cut(phenostoplot$LOD, groups)]

# keep track of the intervals to create a legend
intervals = levels(cut(phenostoplot$LOD, groups, right=T))

# make brackets blank
intervals<-gsub(","," - ",intervals)
intervals<-gsub("\\(","",intervals)
intervals<-gsub("\\]","",intervals)

phenostoplot = modeldf[startsWith(modeldf$phenotype, pattern ="NDVI"),]
check = subset(dfsubset, chromosome == 2)
# reverse order for y axis
# add new column of integers
phenostoplot$ID = seq.int(nrow(phenostoplot))
# sort based on them
phenostoplot = phenostoplot[with(phenostoplot, order(-ID)),]

  filename = addStampToFilename("NDVIstepwiseAllQTLbyDASwithLegend", "pdf")
  pdf(filename, width=17, height=11)
  # fix margins
  # default 
  # par(mar=c(5.1, 4.1, 4.1, 2.1))
  # par(mar=c(1,1,1,1)) # margins: bottom, left, top and right
  plotMultiQTL(cross, 
               stats=phenostoplot, # stats df from stepwise output
               phes=NULL, chrs=NULL, peak=NULL,right=NULL, left=NULL, # if you want to pass individual vectors instead of a combined df
               col=palettevector, # pass a vector to manually color (length of df)
               chr.subset=c(1:12), # subset here
               ylabelcex=1, 
               rugsize=NULL,
               cex=1, 
               pch=19,
               lty=1,
               lwd=2,
               plotQTLdensity=TRUE, 
               binwidth=5, 
               # mark.epi=FALSE, 
               colbychr=F, # set this to F if you want to define custom colors
               palette=rainbow, 
               setmargin=c(8,4,3,2), # adjust legend to fit: bot, left, top, right
               outline=T, background=TRUE
               , main="NDVI stepwise QTL by DAS",
               # xlab= "" # define this afterwards
  )
  # adjust x axis title position
  mtext("Chromosome", side=1, line=2)
  # add legend

  xmin <- par("usr")[1]
  xmax <- par("usr")[2]
  ymin <- par("usr")[3]
  ymax <- par("usr")[4]
  # put it 1/3 of way in middle of plot, and down below the x axis
  legend(xmax/3.3, -(ymax/7), intervals, pch = 19, col = mypalette, title="peak LOD", horiz=T, xpd=NA, cex=0.9)
dev.off()
  

# plot just chr 2 with just NDVI
library(gdata)
phenostoplot = modeldfwithBayes[startsWith(modeldfwithBayes$phenotype, pattern ="NDVI"),]
phenostoplot = subset(phenostoplot, chromosome == 2)
phenostoselect = phenostoplot$phenotype
modelListchr2 = modelList[phenostoselect]

chr2qtl = plot.qtlci(  cross, # cross object
             phenostoplot, # df of stats
             "phenotype", # these are names of cols in df for phenotypes
             "chromosome", # these are names of cols in df for chr
             "lowCIpos", # these are names of cols in df for left CI
             "position", # these are names of cols in df for center pos
             "hiCIpos", # these are names of cols in df for right CI
             models=modelListchr2, # list of models
             plottype="lines", # lines or lodprofile
             standardize=TRUE,
             ci.thickness=1.5,
             point.size=2
)
chr2qtl

filename = addStampToFilename("chr2QTLlinesplot", "pdf")
ggsave(filename, chr2qtl, width=25, height=20, units="in", scale=2, limitsize=F)
