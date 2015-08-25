library("netsim")

ris <- simulatenet(N=5, INdegree="out", times=seq(1,5))
cat("ok")