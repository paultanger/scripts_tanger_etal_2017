# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

# a list of 3 lists for each rep for each pheno
load("finalphenosDS2013_20150518_1230.Robject")

# rename line to ID
colnames(finalphenos)[1] = "ID"

# remove parents (not sure how they got there..)
parents = c("IR64", "Aswina")
finalphenos2 = finalphenos[!(finalphenos$ID %in% parents),]

finalphenos2$pheno_rep = paste(finalphenos2$phenotype, finalphenos2$rep, sep="_")
finalphenos3 = dcast(finalphenos2, ID ~ pheno_rep, value.var="mean", mean, na.rm=T)

# save this for Rqtl in a place we can access via the server
datapath = "~/Desktop/Dropbox/working folder/Rprojects/RILpop/data/"
setwd(datapath)
filename = addStampToFilename("FinalPhenosByRepForRqtl", "csv")
write.table(finalphenos3, file=filename, na="-", col.names=T, row.names=F, sep=",", quote=FALSE)


