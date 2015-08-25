library("netsim")

ris <- simulatenet(N=5, connectivity="random", times=seq(1,5))
cat("ok")