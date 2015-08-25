library(netsim.gui)

p<-read.table("../../parameters.txt")

unlink("moduli.txt")

Dati <- matrix(c(1.793561,1.793561,1.793561,1.793561,1.793561,4.239327,4.239327,4.239327,4.239327,4.239327,1.031338,1.031338,1.031338,1.031338,1.031338,1.789676,1.789676,1.789676,1.789676,1.789676,7.767749,7.767749,7.767749,7.767749,7.767749), 5, 5)

pesi <- matrix(c(0,0,0,-1,0,-1,0,-1,0,0,0,-1,0,0,-1,-1,0,0,0,0,0,0,-1,0,0), 5, 5)

rules <- list(c(5), c(5), c(5), c(5), c(1))

N<-dim(Dati)[1]

ext.W<-matrix(0, N, 2)
ext.W[2,1]<-ext.W[3,1]<-ext.W[2,2]<-1

ris <- simulateprofiles(weights=pesi, EXT.IN=ext.W, EXT.FUN=list("sin(x)", "if(x < 1, 0, x - 1)"), Rules=rules, lambda=p[,1], alpha=p[,2], beta=p[,3], x0=p[,4], Xmin=p[,5], Xmax=p[,6], params=c(0,0,0,0,0,0), times=seq(0, 5, 0.1), save=TRUE)

cat("ok")
