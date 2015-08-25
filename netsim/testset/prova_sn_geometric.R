library("netsim")

ris <- simulatenet(N=5, connectivity="geometric", times=seq(1,5))
cat("ok")