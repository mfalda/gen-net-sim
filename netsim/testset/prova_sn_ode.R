library("netsim")

ris <- simulatenet(N=5, method="rkf45", times=seq(1,5))
cat("ok")