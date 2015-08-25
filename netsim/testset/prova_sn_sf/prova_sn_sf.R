library(netsim)

p<-read.table("../parameters.txt")

init()

simulatenet(N=5, connectivity="scale free", lambda=p[,1], alpha=p[,2], beta=p[,3], x0=p[,4], Xmin=p[,5], Xmax=p[,6], times=seq(0, 5, 0.1))

cat("ok")