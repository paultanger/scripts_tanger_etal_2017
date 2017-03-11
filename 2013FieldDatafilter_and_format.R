# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

finaldata = read.csv("DS2013f7Ril_Masterfile.csv", na.strings="nd",
                     colClasses=c(rep("factor", 3), 
                                  "numeric", "numeric", 
                                  rep("factor", 7), 
                                  rep(NA, 13)))

# remove extra rows and change group to 1-8
finaldata = finaldata[,-c(26:27)]
grouplevels = c(1:8)
levels(finaldata$Group) = grouplevels

# get full list of RILs
AllRILs    = read.csv("SampleListForRwithparents.csv", colClasses=c("factor", "factor", "logical", "logical", "factor", "logical"))

# merge and get more info
colnames(finaldata)[8] = "IRGC_Code"
data = merge(AllRILs, finaldata, by="IRGC_Code", all.y=T)

# get rid of plots with less than 13 plants.. unreliable..
# check how evenly distributed
datatodiscard = data[data$PLTS.PLOT <= 12,]
# write.csv(datatodiscard, "discardedplots.csv")
data = data[data$PLTS.PLOT > 12,]

# separate metadata
as.data.frame(colnames(data))
meta = data[,c(1:17)]

# separate DTF and DTH.. using unique ID to merge later
flowering = data[,c(9,18:19)]

# plant level data.. to get means and merge with ID
plantdata = data[, c(9,20:28)]

# try calculate HI here.. plant level..
plantdata$HI_PLT1 = plantdata$PANWT_PLT1 / (plantdata$FW_PLT1 + plantdata$PANWT_PLT1)
plantdata$HI_PLT2 = plantdata$PANWT_PLT2 / (plantdata$FW_PLT2 + plantdata$PANWT_PLT2)
plantdata$HI_PLT3 = plantdata$PANWT_PLT3 / (plantdata$FW_PLT3 + plantdata$PANWT_PLT3)

# melt data into long format.. this requires reshape, sorry I never updated this code..
detach("package:reshape2", unload=TRUE)
library(reshape)
plantdatalong = melt(plantdata, id.vars=1, variable_name="plantID")

# get rid of plant numbers.. it doesn't matter..
phenos = strsplit(as.character(plantdatalong$plantID), "_")
plantdatalong$pheno = sapply(phenos, "[", 1)

# sort to check that this works
plantdatalong = plantdatalong[with(plantdatalong, order(ID, plantID)), ]

#now need reshape2
detach("package:reshape", unload=TRUE)
library(reshape2)

# get plot means
plotdata.m = dcast(plantdatalong, ID ~ pheno, mean, na.rm=F)
plotdata.m = droplevels(plotdata.m)

# merge plot level data together
DS2013data = merge(meta, flowering, by="ID")
DS2013data = merge(DS2013data, plotdata.m, by="ID")

# per steve's email, fix the row and col.. for some reason they 
# renumber rows and cols for each group, so row and col aren't unique plots
# steve's email:
# The row and column numbers are relative to each Group, with each group having a total of 8 columns and the number of rows depending on the size of the group.  The table also has a Field number ending in .1 or .2 for each side of the field.  So if you prefer to have row/column numbers based on the entire field you could simply order the data by Field,Group, Row and create a new column with row numbers extending the length of the field and modify the column numbers for plots on one side of each field (ie., 519.2) to be 9-16. 

# sort by field, group, row
DS2013data = DS2013data[with(DS2013data, order(Field, Group, Row)), ]

# just keep stuff we need.. we can always get other stuff later..
DS2013data = DS2013data[, c(1,3,9,10,11,14,16,18:23)]

# rename cols
colnames(DS2013data) = c("ID", "line", "field", "row", "col", "rep", "group", "DTF", "DTH", "biomass", "HI", "height", "grain")

# drop NA
DS2013data = DS2013data[!is.na(DS2013data$line),]

# add a number to the rows for the right groups, depending on the rep
# this is kind of a nightmare because different reps have different rows 
# in each group and need to refer to the field layout to determine
# the correct groups to modify

# create a new row to modify
DS2013data$fixedrow = DS2013data$row

# get rid of factors
DS2013data$group = as.numeric(as.character(DS2013data$group))
DS2013data$rep = as.numeric(as.character(DS2013data$rep))

# rep 1
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 8 & DS2013data$rep == 1, # based on these conditions
         DS2013data$fixedrow + 25, # change it 
         DS2013data$fixedrow) # otherwise do nothing

DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 5 & DS2013data$rep == 1, # based on these conditions
         DS2013data$fixedrow + 24, # change it 
         DS2013data$fixedrow) # otherwise do nothing
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 3 & DS2013data$rep == 1, # based on these conditions
         DS2013data$fixedrow + 24, # change it 
         DS2013data$fixedrow) # otherwise do nothing
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 2 & DS2013data$rep == 1, # based on these conditions
         DS2013data$fixedrow + 25, # change it 
         DS2013data$fixedrow) # otherwise do nothing
# rep 2
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 5 & DS2013data$rep == 2, # based on these conditions
         DS2013data$fixedrow + 24, # change it 
         DS2013data$fixedrow) # otherwise do nothing
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 2 & DS2013data$rep == 2, # based on these conditions
         DS2013data$fixedrow + 24, # change it 
         DS2013data$fixedrow) # otherwise do nothing
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 4 & DS2013data$rep == 2, # based on these conditions
         DS2013data$fixedrow + 25, # change it 
         DS2013data$fixedrow) # otherwise do nothing
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 7 & DS2013data$rep == 2, # based on these conditions
         DS2013data$fixedrow + 25, # change it 
         DS2013data$fixedrow) # otherwise do nothing
# rep 3
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 4 & DS2013data$rep == 3, # based on these conditions
         DS2013data$fixedrow + 24, # change it 
         DS2013data$fixedrow) # otherwise do nothing
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 3 & DS2013data$rep == 3, # based on these conditions
         DS2013data$fixedrow + 24, # change it 
         DS2013data$fixedrow) # otherwise do nothing
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 8 & DS2013data$rep == 3, # based on these conditions
         DS2013data$fixedrow + 25, # change it 
         DS2013data$fixedrow) # otherwise do nothing
DS2013data$fixedrow = # reassign col
  ifelse(DS2013data$group == 2 & DS2013data$rep == 3, # based on these conditions
         DS2013data$fixedrow + 25, # change it 
         DS2013data$fixedrow) # otherwise do nothing

# now fix columns.. add 8 to correct groups
rep1groups = c(7,5,4,2)
rep2groups = c(8,2,3,7)
rep3groups = c(1,3,6,2)

DS2013data$fixedcol = DS2013data$col

DS2013data$fixedcol = # reassign col
  ifelse(DS2013data$group %in% rep1groups & DS2013data$rep == 1, # based on these conditions
         DS2013data$fixedcol + 8, # change it 
         DS2013data$fixedcol) # otherwise do nothing

DS2013data$fixedcol = # reassign col
  ifelse(DS2013data$group %in% rep2groups & DS2013data$rep == 2, # based on these conditions
         DS2013data$fixedcol + 8, # change it 
         DS2013data$fixedcol) # otherwise do nothing

DS2013data$fixedcol = # reassign col
  ifelse(DS2013data$group %in% rep3groups & DS2013data$rep == 3, # based on these conditions
         DS2013data$fixedcol + 8, # change it 
         DS2013data$fixedcol) # otherwise do nothing

# as csv.. with dots for NA
# also, fields have to be one field (not .1 and .2)
DS2013data$field = as.integer(DS2013data$field)

# # sort by line and rep
DS2013data = DS2013data[with(DS2013data, order(line, rep)), ]
# remove old row col
DS2013data = DS2013data[,-c(4,5)]
DS2013data = DS2013data[,c(1:5,12:13,6:11)]

colnames(DS2013data)[colnames(DS2013data) == "fixedrow"] = "row"
colnames(DS2013data)[colnames(DS2013data) == "fixedcol"] = "col"

fileprefix = paste("DS2013phenofinalFIXROWCOLfixfield_", format(Sys.time(),"%Y%m%d_%H%M"), sep="")
# write.table(DS2013data, file=paste(fileprefix, "tsv", sep="."),row.names=F, sep="\t", quote=F)
