library("netsim")

ris <- simulatenet(N=5, act.fun="linear", times=seq(1,5))
cat("ok")