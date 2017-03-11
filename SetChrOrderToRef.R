# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

SetChrOrderToRef = function(cross){
  
  tofix = vector()
  # for each chr
  for (i in names(cross$geno)){
    print (i)
    print(names(cross$geno[[i]]$map))
    # get bp of first and last marker in map
    first = sapply(strsplit( names(cross$geno[[i]]$map[1]), "_"), "[[", 2) 
    last = sapply(strsplit( names(cross$geno[[i]]$map[length(names(cross$geno[[i]]$map))]), "_"), "[[", 2) 
    # change to numbers
    first = as.numeric(first)
    last = as.numeric(last)
    print(first)
    print(last)
    # if last is smaller than first, it needs to be reversed
    if (last < first) {
      # add to vector
      tofix = c(tofix, i)    
    }
  }
  
  # fix them
  cross2 = cross
  for (i in tofix){
    print(i)
    cross2 = switch.order(cross2, chr=i, nmar(cross)[i]:1)
  }
  
  # return new cross object
  return(cross2)
}