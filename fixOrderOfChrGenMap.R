# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

setwd("../input")

load("crossallphenos_20150521_1444.Robject")

cross2 = SetChrOrderToRef(cross)

# check that it worked
crosscheck = SetChrOrderToRef(cross2)

# we need to re-run calc geno
cross = calc.genoprob(cross2, map.function="kosambi", step=1)

# save it
filename = addStampToFilename("crossallphenosOrderFixedReRunCalcGeno", "Robject")
save(cross, file=filename)
