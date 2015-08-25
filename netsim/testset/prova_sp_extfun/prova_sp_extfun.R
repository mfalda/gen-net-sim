library(netsim)

p<-read.table("../parameters.txt")

init()

Dati <- matrix(c(1.793561,1.793561,1.793561,1.793561,1.793561,4.239327,4.239327,4.239327,4.239327,4.239327,1.031338,1.031338,1.031338,1.031338,1.031338,1.789676,1.789676,1.789676,1.789676,1.789676,7.767749,7.767749,7.767749,7.767749,7.767749), 5, 5)

pesi <- matrix(c(0,0,0,-1,0,-1,0,-1,0,0,0,-1,0,0,-1,-1,0,0,0,0,0,0,-1,0,0), 5, 5)

rules <- list(c(5), c(5), c(5), c(5), c(1))

N<-dim(Dati)[1]

ext.W<-matrix(0, N, 2)
ext.W[2,1]<-ext.W[3,1]<-ext.W[2,2]<-1

f1<-function(t){y<-sin(t); return(y)}
f2<-function(t){if (t<1) y<-0 else y<-t-1; return(y)}

ris <- simulateprofiles(weights=pesi, EXT.IN=ext.W, EXT.FUN=list(f1, f2), Rules=rules, lambda=p[,1], alpha=p[,2], beta=p[,3], x0=p[,4], Xmin=p[,5], Xmax=p[,6], times=seq(0, 5, 0.1))

write.table(ris, file="SIMdata_sp.txt")

cat("ok")
