# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

ExtractGenoFromCross = function(qtlobject) {
  # extract genotype data from object to use for joinmap, MSTmap, or our CS tool
  genos = as.data.frame(qtlobject$geno[[1]]$data)
  i = 1
  for (i in 2:length(qtlobject$geno)) {
    genos = cbind(as.matrix(genos), qtlobject$geno[[i]]$data) 
  }
    
  # what about pull.geno? or write.cross
  
  # change back to A and B
  genos[genos == 1] = "A"
  genos[genos == 2] = "B"
  genos[is.na(genos)] = "-"
  
  genos = as.data.frame(genos)
  
  # add line names back
  finalgenos = cbind(qtlobject$pheno[1], genos)
  
  geno_joinmap = finalgenos
  
  # transpose
  geno.t = t(geno_joinmap)
  
  # add col names
  colnames(geno.t) <- geno_joinmap[,1]
  
  # remove first row
  geno.t = geno.t[-1,]
  
  # make row.names into column
  geno.t2 = as.data.frame(cbind(line = rownames(geno.t), geno.t)) 
  
  # get map
  themap = pull.map(qtlobject, as.table=T)
  
  # merge map distances and geno data
  genomap = merge(themap, geno.t2, by="row.names")
  
  # sort
  genomap2 = with(genomap, genomap[order(chr, pos), ])
  
  # transpose
  genomap = t(genomap2)
  
  # add line col as "ID"
  genomap <- cbind(ID = rownames(genomap), genomap)
  
  # rename cols
  colnames(genomap) <- as.character(genomap[1,])
  
  # fix ID col
  colnames(genomap)[1] = "ID"
  # remove old row with col names
  genomap = genomap[-c(1,4),]
  
  # remove chrom and marker and cM from ID col
  genomap[1,1] = ""
  genomap[2,1] = ""
  
  return(genomap)
}