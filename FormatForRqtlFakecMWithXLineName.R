formatforRqtlfakecM = function(genotypes){
  # remove some cols
  genotypes = genotypes[,-c(2)]
  
  # create column with fake cM
  genotypes$pos = as.numeric(as.character(genotypes$pos)) / 10000
  # transpose
  genotypes = t(genotypes)
  
  # add line col as "ID"
  genotypes <- cbind(ID = rownames(genotypes), genotypes)
  
  # remove chrom and marker and cM from ID col
  genotypes[1,1] = ""
  genotypes[2,1] = ""
  genotypes[3,1] = ""

  # rename cols
  colnames(genotypes) <- as.character(genotypes[1,])
  # fix ID col
  colnames(genotypes)[1] = "ID"
  # remove old row with col names
  genotypes = genotypes[-1,]

  # change line names to match phenotype data (just number)
  # split = strsplit(genotypes[,1], "_")
  # split = sapply(split, "[", 2)
  # use character vector as new col names
  
  # remove X from line names
  split = strsplit(genotypes[,1], "X")
  split = sapply(split, "[", 2)
  genotypes = cbind(as.data.frame(as.character(split)), genotypes)
  genotypes = cbind(split, genotypes)
  
  genotypes = genotypes[,-2]
  colnames(genotypes)[colnames(genotypes)=="split"] = "ID"
  

  # genotypes$ID = split
  genotypes[1,1] = ""
  genotypes[2,1] = ""
  
  # try with numeric instead of character.. na as ""
  genotypes$ID = as.numeric(genotypes$ID)
  
  return(genotypes)
}