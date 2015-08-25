library("netsim")

ris <- simulatenet(N=50, save=TRUE, itera=5, times=seq(1,10))
cat("ok")