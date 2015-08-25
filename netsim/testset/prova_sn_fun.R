library("netsim")

ris <- simulatenet(N=5, f.pr.and="1 / (1 + e ^ (-10 * (x - 0.5)))", times=seq(1,5))
cat("ok")