# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

FormatForJoinMap = function(qtlobject) {
  # extract genotype data from object to use for joinmap, MSTmap, or our CS tool
  genos = as.data.frame(qtlobject$geno[[1]]$data)
  i = 1
  for (i in 2:length(qtlobject$geno)) {
    genos = cbind(as.matrix(genos), qtlobject$geno[[i]]$data) 
  }
    
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
  
  # add dumb column for joinmap
  geno_joinmap = cbind(marker = geno.t2[,1], class="(a,b)", geno.t2[,2:ncol(geno.t2)])
  
  # check stuff
  individuals = ncol(geno_joinmap)-2
  individuals
  loci = nrow(geno_joinmap)
  loci
  
  return(geno_joinmap)
}