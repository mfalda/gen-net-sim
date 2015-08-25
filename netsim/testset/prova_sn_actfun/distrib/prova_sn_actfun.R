library(netsim.gui)

p<-read.table("../../parameters.txt")

unlink("moduli.txt")

simulatenet(N=5, act.fun="linear", lambda=p[,1], alpha=p[,2], beta=p[,3], x0=p[,4], Xmin=p[,5], Xmax=p[,6], params=c(0,0,0,0,0,0), times=seq(0, 5, 0.1), save=TRUE)

cat("ok")