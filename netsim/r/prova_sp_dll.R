setwd("D:\\Progetti R\\netsim_hg\\sim")

dyn.load("simulatenet.dll")
dyn.load("simulateprofiles.dll")

source("simulatenet.R")
source("simulateprofiles.R")

pesi <- matrix(c(0,0,0,-1,0,-1,0,-1,0,0,0,-1,0,0,-1,-1,0,0,0,0,0,0,-1,0,0), 5, 5)

lambda <-c(1.1056752,0.9747519,1.0007163,1.0061211,1.1354204)

alpha <- c(9.869912,9.955962,10.056057,9.922354,9.950865)

rules <- list(c(5), c(5), c(5), c(5), c(1))

beta <- c(0.4934634,0.4930697,0.4926102,0.4976155,0.4905647)

times<-seq(1,5)
save <- FALSE
itera <- 1

ris <- simulateprofiles(weights=pesi, Rules=rules, lambda=lambda, alpha=alpha, method="rkf45", beta=beta, times=times, save=F)

cat("ok")

net<-simulatenet(N=5,connectivity="MTM",act.fun="sigmoidal", method="rkf45", save=F)
W.matrix<-net[[1]][[2]]
R<-net[[2]]
param.lambda<-net[[1]][[3]]
param.alpha<-net[[1]][[4]]
param.beta<-net[[1]][[5]]
baseline<-runif(5,0,1)
f1<-"sin(t)"
f2<-"if(t<1, 0, -t-1)"
ext.W<-matrix(0,5,2)
ext.W[2,1]<-ext.W[3,1]<-ext.W[2,2]<-1
Data.new<-simulateprofiles(X0=baseline,Xmax=rep(10,5),lambda=param.lambda,weights=W.matrix,Rules=R,act.fun="sigmoidal",alpha=param.alpha, beta=param.beta, method="rkf45", EXT.IN=ext.W, EXT.FUN=list(f1, f2), save=F)

cat("ok")
