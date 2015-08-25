args = commandArgs(TRUE)
N = 5

ambiente <- new.env()
assign("param_writetable", 0)
assign("param_checkconn", 1)
assign("param_createNEG", 1)
assign("param_createRules", 1)
assign("param_target", 1)
assign("param_connectivitygeometric", 1)
assign("param_connectivityscalefree", 1)
	assign("param_Score", 1)
	assign("param_Scoresf", 1)
assign("param_connectivitymodular", 1) # dipende anche da score
	assign("param_module1", 1)
	assign("param_module2", 1)
	assign("param_module3", 1)
	assign("param_assignnodes", 1)
	assign("param_clustercoeff", 1)
assign("param_connectivityrandom", 1)

library(netsim1all)
#~ f1<-function(x) {y<-1/(1+exp(-10*(x-0.5))); return(y)}
net<-simulatenet(N,connectivity="MTM",f.pr.and="1 / (1 + e ^ (-10 * (x - 0.5)))", act.fun="sigmoidal", method="lsoda")


Data<-net[[1]]
W.matrix<-net[[2]]
R<-net[[3]]
param<-net[[4]]
N<-dim(Data)[1]
baseline<-runif(N,0,1)
#~ f1<-function(t){y<-sin(t); return(y)}
#~ f2<-function(t){if (t<1) y<-0 else y<-t-1; return(y)}
ext.W<-matrix(0, N, 2)
ext.W[2,1]<-ext.W[3,1]<-ext.W[2,2]<-1
Data.new<-simulateprofiles(x0=baseline,Xmax=rep(10,N),lambda=param[,1],weights=W.matrix,Rules=R,act.fun="sigmoidal",alpha=param[,2], beta=param[,3], method="lsoda", EXT.IN=ext.W, EXT.FUN=list("sin(x)", "if(x < 1, 0, x - 1)"))

