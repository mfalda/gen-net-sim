args = commandArgs(TRUE)
N = as.integer(args[1])
N <- 200

library(netsim)
net<-simulatenet(N, connectivity="MTM", method="rkf45")

Data<-net[[1]][[1]]
W.matrix<-net[[1]][[2]]
R<-net[[2]]
lambda<-net[[1]][[3]]
alpha<-net[[1]][[4]]
beta<-net[[1]][[5]]
N<-dim(Data)[1]
baseline<-runif(N,0,1)
Data.new<-simulateprofiles(x0=baseline,Xmax=rep(10,N),lambda=lambda,weights=W.matrix,Rules=R,act.fun="sigmoidal",alpha=alpha, beta=beta, method="rkf45")
