# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.


mergeparentalcalls = function (SNPcalls) {
  
  # Aswina
  for (i in 1:nrow(SNPcalls)){
    # if all three equal
    if (SNPcalls[i,5]==SNPcalls[i,6] && SNPcalls[i,6]==SNPcalls[i,7]){
      SNPcalls$Aswinam[i] <- SNPcalls[i,5]
    }
    # if first 2 equal
    else if (SNPcalls[i,5]==SNPcalls[i,6]){
      SNPcalls$Aswinam[i] <- SNPcalls[i,5]
    }  
    # if last 2 equal
    else if (SNPcalls[i,6]==SNPcalls[i,7]){
      SNPcalls$Aswinam[i] <- SNPcalls[i,6] 
    }
    # if first & last equal
    else if (SNPcalls[i,5]==SNPcalls[i,7]){
      SNPcalls$Aswinam[i] <- SNPcalls[i,5]
    }
  }
  
  # check it
  #test <- hap[, c(5:7, which(colnames(hap)=="Aswina"))]
  #test
  
  # IR64
  for (i in 1:nrow(SNPcalls)){
    # if all three equal
    if (SNPcalls[i,8]==SNPcalls[i,9] && SNPcalls[i,9]==SNPcalls[i,10]){
      SNPcalls$IR64m[i] <- SNPcalls[i,8]}
    # if first 2 equal
    else if (SNPcalls[i,8]==SNPcalls[i,9]){
      SNPcalls$IR64m[i] <- SNPcalls[i,8]
    }  
    # if last 2 equal
    else if (SNPcalls[i,9]==SNPcalls[i,10]){
      SNPcalls$IR64m[i] <- SNPcalls[i,9] 
    }
    # if first & last equal
    else if (SNPcalls[i,8]==SNPcalls[i,10]){
      SNPcalls$IR64m[i] <- SNPcalls[i,8]
    }
  }
  
  return (SNPcalls)
}