# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

load("DS2013AllMeansNewWithCTDFixHI.Robject")

DS2013AllMeans$phenotype = as.factor(DS2013AllMeans$phenotype)
levels(DS2013AllMeans$phenotype)

# change order
DS2013AllMeans$phenotype  <- factor(DS2013AllMeans$phenotype , levels = 
                                      c("DTF",
                                        "DTH",
                                        "biomass",
                                        "height",
                                        "grain",
                                        "HI",
                                        "CTD",
                                        "Chla",
                                        "HTPheight",
                                        "NDRE",
                                        "NDVI",
                                        "NIR_760",
                                        "R_670",
                                        "RE_730"))

# rename grain factor
levels(DS2013AllMeans$phenotype) <- sub("^grain$", "grain_yield", levels(DS2013AllMeans$phenotype))

# define final phenos to include with all
finalphenos = c("DTH", "biomass", "height", "grain_yield", "HI")

CanTemp   = subset(DS2013AllMeans, phenotype %in% c("CTD", finalphenos))
Chla      = subset(DS2013AllMeans, phenotype %in% c("Chla", finalphenos))
HTPheight = subset(DS2013AllMeans, phenotype %in% c("HTPheight", finalphenos))
NDRE      = subset(DS2013AllMeans, phenotype %in% c("NDRE", finalphenos))
NDVI      = subset(DS2013AllMeans, phenotype %in% c("NDVI", finalphenos))

finalstuff = subset(DS2013AllMeans, phenotype %in% c(finalphenos) & date == c("final") )

# skip all this stuff below.. just make a list of 1..
finalphenoslist = list(finalstuff)

# make a list
phenotypes = list(CanTemp, Chla, HTPheight, NDRE, NDVI)

# name items in list
phenoNames = c("CTD", "Chla", "HTPheight", "NDRE", "NDVI")
names(phenotypes) = phenoNames

# drop levels
phenotypes = lapply(phenotypes, function(x) droplevels(x))
library(reshape2)
phenotypesforCorr = lapply(phenotypes, function(x) {
  dcast(x, line + rep ~ phenotype + dateformat, value.var="mean", mean, na.rm=T)
})

# prep for correlations
library(reshape2)

# set colors with alpha
# colors lovell used for other figs
cols<- c("#000000", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
# convert these:"#56B4E9", "#009E73", "#F0E442"
alpha = 0.3
# have to divide by 255...
# rep 1 blue
# rep 2 red
# rep 3 yellow
alphablue = rgb(70/255, 172/255, 233/255, alpha, "alphablue")
alphared = rgb(180/255, 16/255, 32/255, alpha, "alphared")
alphagreen = rgb(226/255, 210/255, 2/255, alpha, "alphagreen")
colors = c(alphablue, alphared, alphagreen)

# define function to put the r values in the top of the scatterplot
#http://gettinggeneticsdone.blogspot.com/2011/07/scatterplot-matrices-in-r.html
# install.packages("Hmisc")
library(Hmisc)
panel.cor <- function(x, y, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x,y, use="pairwise.complete.obs") # na.or.complete doesn't work correctly.. removes rows..
  p <- rcorr(x,y, type="pearson")$P[1,2]
  rval <- formatC(r, 2, format="f")
  rval <- paste("r = ", rval, sep="")
  pval <- formatC(p, 2, format="f")
  pval <- paste("p = ", pval, sep="")
  if(pval == "p = 0.00"){ pval = "p < 0.01"}
  text(0.5, 0.7, rval, cex=1.5) # if this is too big, it messes up pdf export
  text(0.5, 0.3, pval, cex=1.5)
}

# filename = addStampToFilename("scatterplots", "pdf")
setwd(datapath)
i= 1
# #   postscript(filename, width=10, height=8.5, paper="special", family="Helvetica", fonts="Helvetica", horizontal=F, onefile=F, bg="white") # haven't tried this..
# # pdf(filename, width=15, height=15) # file is huge!
# for(i in i:length(phenotypes)){
#     # this seems to work the best
#     filename = addStampToFilename(paste0(names(phenotypesforCorr[i]), "_scatterplot"), "jpg")
#     jpeg(filename, width=15, height=15, units="in", res=150, quality=90)
#     # eps instead
#     # filename = addStampToFilename(paste0(names(phenotypesforCorr[i]), "_scatterplot"), "eps")
#     # postscript(filename, width=15, height=15, paper="special", family="Helvetica", fonts="Helvetica", horizontal=F, onefile=F, bg="white")
#     # make the plot
#   pairs(phenotypesforCorr[[i]][c(3:length(colnames(phenotypesforCorr[[i]])))],
#      #   main=paste0("line means of ", names(phenotypesforCorr[i])," , reps: 1=red, 2=green, 3=blue"),
#   #       col="#00000033",
#   #       pch=21,
#         col = c(colors)[unclass(phenotypesforCorr[[i]]$rep)],
#         cex = 0.5
#         ,upper.panel=panel.cor
#         )
#   dev.off()
# }


# or plot separately for each rep...
# can't do rep3.. have to do separately below..
i= 1
for(i in i:length(phenotypes)){
  for(j in levels(phenotypesforCorr[[i]]$rep)){
    # name the file
    filename = addStampToFilename(paste0(names(phenotypesforCorr[i]), "_rep_", j, "_scatterplot"), "jpg")
    jpeg(filename, width=15, height=15, units="in", res=150, quality=90)
    # subset data
    tempdata = subset(phenotypesforCorr[[i]], rep == as.character(j))
    # make the plot
    # skip grain and 4/17 data since missing some reps
    pairs(tempdata[c(3:5,8:(length(colnames(tempdata))-1))],
    # pairs(tempdata[c(3:6,8:(length(colnames(tempdata))))],
          main=paste0("line means of ", names(phenotypesforCorr[i])," , Cohort ", j),
          cex = 0.5, 
          upper.panel=panel.cor
    )
    dev.off()
  }
}


i= 1
for(i in i:length(phenotypes)){
  # name the file
  filename = addStampToFilename(paste0(names(phenotypesforCorr[i]), "_rep_", 3, "_scatterplot"), "jpg")
  jpeg(filename, width=15, height=15, units="in", res=150, quality=90)
  # subset data
  tempdata = subset(phenotypesforCorr[[i]], rep == 3)
  # make the plot
  pairs(tempdata[c(3:7,8:(length(colnames(tempdata))))],
        main=paste0("line means of ", names(phenotypesforCorr[i])," , Cohort ", 3),
        cex = 0.5, 
        upper.panel=panel.cor
  )
  dev.off()
}