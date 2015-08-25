library("netsim")

ris <- simulatenet(N=5, connectivity="scale free", times=seq(1,5))
cat("ok")