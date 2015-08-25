library("netsim")

ris <- simulatenet(N=5, connectivity="MTM", times=seq(1,5))
cat("ok")