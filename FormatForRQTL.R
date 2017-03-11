# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

FormatForRqtl = function(data){
  
  # transpose it
  datat = t(data)
#   head(geno_rqtl.t[1:5,1:5])
  
  # add line col as "ID"
  datat <- cbind(ID = rownames(datat), datat)
  
  # remove chrom and marker and cM from ID col
  datat[1,1] = ""
  datat[2,1] = ""
  datat[3,1] = ""
#   datat[1:5,1:5]
  
  # rename cols
  colnames(datat) <- as.character(datat[1,])
  # fix ID col
  colnames(datat)[1] = "ID"
  # remove old row with col names
  datat = datat[-1,]
#   str(geno_rqtl.t)
  
  # remove X from IDs
#   split = strsplit(datat[,1], "X")
#   split = sapply(split, "[", 2)
#   # use character vector as new col names
#   datat = cbind(as.data.frame(as.character(split)), datat)
#   datat = datat[,-2]
#   colnames(datat)[colnames(datat)=="as.character(split)"] = "ID"
#   datat$ID = as.character(datat$ID)
#   datat[1,1] = ""
#   datat[2,1] = ""
  
  # try with numeric instead of character.. na as ""
#   datat$ID = as.numeric(datat$ID)
  
  return (datat)
}