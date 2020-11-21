library(propr)
library(zCompositions)
data(caneToad.counts)

counts = caneToad.counts[,1:10000]
counts=cmultRepl(counts,method="CZM",label=0,output="p-counts")
counts2 = caneToad.counts[,1:1000]

clr = propr:::proprCLR(counts)
clr2 = propr:::proprCLR(counts2)

pro = propr(clr, "rho", p=20, ivar=NA)
pro2 = propr(clr2, "rho", p=20, ivar=NA)