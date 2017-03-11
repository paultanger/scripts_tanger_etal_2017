# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

# reports logical cores (generally 2x physical cores)
# detectCores()
# use n - 1 cores
numworkers = detectCores() - 1
message("using up to ", numworkers, " cores")

runDropOneMarker = function(crossobject, chrs, maptype){
  result = mclapply(chrs, droponemarker, cross=crossobject, map.function=maptype, mc.cores=numworkers)
  # this is a list, so name the items and combine them
#   names(result) = chrs
  df <- do.call("rbind", result)
  return(df)
}