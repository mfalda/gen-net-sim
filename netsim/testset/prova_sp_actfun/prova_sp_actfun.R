library(netsim)

p<-read.table("../parameters.txt")

init()

Dati <- matrix(c(1.793561,1.793561,1.793561,1.793561,1.793561,4.239327,4.239327,4.239327,4.239327,4.239327,1.031338,1.031338,1.031338,1.031338,1.031338,1.789676,1.789676,1.789676,1.789676,1.789676,7.767749,7.767749,7.767749,7.767749,7.767749), 5, 5)

pesi <- matrix(c(0,0,0,-1,0,-1,0,-1,0,0,0,-1,0,0,-1,-1,0,0,0,0,0,0,-1,0,0), 5, 5)

rules <- list(c(5), c(5), c(5), c(5), c(1))

N<-dim(Dati)[1]

ris <- simulateprofiles(weights=pesi, act.fun="linear", Rules=rules, lambda=p[,1], alpha=p[,2], beta=p[,3], x0=p[,4], Xmin=p[,5], Xmax=p[,6], times=seq(0, 5, 0.1))

write.table(ris, file="SIMdata_sp.txt")

cat("ok")
