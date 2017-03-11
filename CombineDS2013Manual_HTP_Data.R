# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

DS2013data = read.delim("DS2013phenofinalFIXROWCOLfixfield_20160929_1943.tsv")

setwd("CS_Final/DONOTTOUCH_NEW/")
file_list = list.files()
# combine them
CSfiles <- do.call("rbind",lapply(file_list,
                                  FUN=function(files){
                                    # get data
                                    df = read.table(files, header=TRUE, sep=",")
                                    # add date column
                                    df$date = substr(files, 3, nchar(files)-7)
                                    # return it to rbind
                                    return(df)
                                  }))

# make factors
CSfiles$ID_Plot = as.factor(CSfiles$ID_Plot) 
CSfiles$date = as.factor(CSfiles$date) 

# do the same for CC files
setwd("../../CC_Final/DONOTTOUCH/")
file_list = list.files()
# combine them
CCfiles <- do.call("rbind",lapply(file_list,
                                  FUN=function(files){
                                    # get datafile
                                    df = read.table(files, header=TRUE, sep=",", stringsAsFactors=F)
                                    # add date column
                                    df$newdate = substr(files, 3, nchar(files)-7)
                                    # return it to rbind
                                    return(df)
                                  }))

# make factors
CCfiles$ID_Plot = as.factor(CCfiles$ID_Plot) 
CCfiles$date = as.factor(CCfiles$newdate) 

# clean up
as.data.frame(colnames(CSfiles))
CSfiles = CSfiles[,c(16,19,3,14,20)]
colnames(CSfiles) = c("rep", "CTD", "ID_plot", "HTPheight", "date")

CCfiles = CCfiles[,c(10,3,13:15,17,16,18,21)]
colnames(CCfiles)[1] = "rep"

# convert to cm
# rep as factor
CSfiles$rep = as.factor(CSfiles$rep) 
CCfiles$rep = as.factor(CCfiles$rep) 

# this is raw data of the sensor over each plot.. 
HTPheightdata = CSfiles
HTPspectraldata = CCfiles

# need to make this long first
HTPheightTemplong = melt(HTPheightdata, id=c("ID_plot", "rep", "date"), na.rm=T, value.name="pheno")

HTPheightTempsummary <- aggregate(pheno ~ ID_plot + rep + variable + date, data=HTPheightTemplong, FUN=function(x) c(length(x), mean(x), sd(x), sd(x)/sqrt(length(x)), sd(x)/mean(x)))
HTPheightTempsummary <- cbind(HTPheightTempsummary[,1:4], as.data.frame(HTPheightTempsummary[,5]))
names(HTPheightTempsummary) <- c("plot", "rep", "pheno", "date", "n", "mean", "SD", "SE", "CV")

# need to make this long first
HTPspectraldatalong = melt(HTPspectraldata, id=c("ID_Plot", "rep", "date"), na.rm=T, value.name="spectra")

HTPspectralsummary <- aggregate(spectra ~ ID_Plot + rep + variable + date, data=HTPspectraldatalong, FUN=function(x) c(length(x), mean(x), sd(x), sd(x)/sqrt(length(x)), sd(x)/mean(x)))
HTPspectralsummary <- cbind(HTPspectralsummary[,1:4], as.data.frame(HTPspectralsummary[,5]))
names(HTPspectralsummary) <- c("plot", "rep", "spectra", "date", "n", "mean", "SD", "SE", "CV")

# add rep date spectra factor
HTPspectralsummary$date_rep_spectra = paste(HTPspectralsummary$date, HTPspectralsummary$rep, HTPspectralsummary$spectra, sep="_")
HTPheightTempsummary$date_rep_pheno = paste(HTPheightTempsummary$date, HTPheightTempsummary$rep, HTPheightTempsummary$pheno, sep="_")

# discard data where less than 4 readings per plot
HTPspectralsummary2 = HTPspectralsummary[HTPspectralsummary$n > 3, ]
HTPheightsummary2 = HTPheightTempsummary[HTPheightTempsummary$n > 3, ]

# discard data where CV > 50 %
HTPspectralsummary3 = HTPspectralsummary2[HTPspectralsummary2$CV < 0.5, ]
HTPheightsummary3 = HTPheightsummary2[HTPheightsummary2$CV < 0.5, ]

# for now, try getting means for each line / date / rep
plotToLine = DS2013data[,c(1:2)]

# add line info
HTPspectralsummary4 = merge(plotToLine, HTPspectralsummary3[,c(1:4,6)], by.x="ID", by.y="plot") # include all.y=T to run next step
HTPheightsummary4 = merge(plotToLine, HTPheightsummary3[,c(1:4,6:8)], by.x="ID", by.y="plot") # include all.y=T to run next step

# now we have the line, so get the means for each date and spectra..
# this shouldn't be that different since we don't have many duplicate plots
HTPspectralsummarylinerepdate <- aggregate(mean ~ line + rep + date + spectra, data=HTPspectralsummary4, FUN=function(x) c(length(x), mean(x), sd(x), sd(x)/sqrt(length(x)), sd(x)/mean(x)))
HTPspectralsummarylinerepdate <- cbind(HTPspectralsummarylinerepdate[,1:4], as.data.frame(HTPspectralsummarylinerepdate[,5]))
names(HTPspectralsummarylinerepdate) <- c("line", "rep", "date", "spectra", "n", "mean", "SD", "SE", "CV")

HTPheightsummarylinerepdate <- aggregate(mean ~ line + rep + pheno + date, data=HTPheightsummary4, FUN=function(x) c(length(x), mean(x), min(x), max(x), sd(x), sd(x)/sqrt(length(x)), sd(x)/mean(x)))
HTPheightsummarylinerepdate <- cbind(HTPheightsummarylinerepdate[,1:4], as.data.frame(HTPheightsummarylinerepdate[,5]))
names(HTPheightsummarylinerepdate) <- c("line", "rep", "pheno", "date", "n", "mean", "min", "max", "SD", "SE", "CV")

# split into canTemp and HTPheight
HTPheight = subset(HTPheightsummarylinerepdate, pheno == "HTPheight")
HTPcanTemp = subset(HTPheightsummarylinerepdate, pheno == "CTD")

# add ID col
HTPheight$spectra = "HTPheight"
HTPcanTemp$spectra = "CTD"

# keep stuff we need
HTPcanTemp = HTPcanTemp[,c(1,2,4,12,6)]
HTPheight = HTPheight[,c(1,2,4,12,6:8)]

# combine spectra and cantemp
HTPmeans = rbind(HTPcanTemp, HTPspectralsummarylinerepdate[,c(1:4,6)])
HTPmeans$spectra = as.factor(HTPmeans$spectra)

# for now, make 3 files.. only thing that is diff is whether using min max or mean of height
HTPmeansHeightMean = rbind(HTPheight[,c(1:4,5)], HTPmeans)

# load manual data

colnames(HTPmeansHeightMean)[4] = "phenotype"

DS2013data$date = as.factor("final")
# make it long
DS2013datal = melt(DS2013data[,c(2,4,8:14)], id=c("line", "rep", "date"), na.rm=T, variable.name="phenotype", value.name="mean")
DS2013AllMeans = rbind(HTPmeansHeightMean, DS2013datal)

# add rep date spectra
DS2013AllMeans$date_rep_phenotype = paste(DS2013AllMeans$date, DS2013AllMeans$rep, DS2013AllMeans$phenotype, sep="_")

# in order to look at correlations between dates for each phenotype, make the data wide

# change date to date
DS2013AllMeans$dateformat = as.Date(DS2013AllMeans$date, format = "%m%d%y")
# format nicely..
DS2013AllMeans$dateformat = format(DS2013AllMeans$dateformat, format= "%m%d%y")
# sort
DS2013AllMeans = DS2013AllMeans[with(DS2013AllMeans, order(dateformat)), ]

AllmeansWide = dcast(DS2013AllMeans, line + rep ~ phenotype + dateformat, value.var="mean", mean, na.rm=T)

filename = addStampToFilename("AllmeansForCorr", "tsv")
# write.table(AllmeansWide, file=filename, sep="\t", quote=F, row.names=F)

setwd("../../")
save(DS2013AllMeans, file="DS2013AllMeansNewWithCTDFixHI.Robject")
save(AllmeansWide, file="AllmeansWideNew.Robject")
