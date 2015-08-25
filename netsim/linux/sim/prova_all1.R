args = commandArgs(TRUE)
N = as.integer(args[1])

ambiente <- new.env()
assign("param_writetable", 1)
assign("param_dinamica", 1)
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
	assign("param_clustercoeff",1)
assign("param_connectivityrandom", 1)

library(netsim1all)
net<-simulatenet(N, connectivity="MTM")

Data<-net[[1]]
W.matrix<-net[[2]]
R<-net[[3]]
param<-net[[4]]
N<-dim(Data)[1]
baseline<-runif(N,0,1)
Data.new<-simulateprofiles(x0=baseline,Xmax=rep(10,N),lambda=param[,1],weights=W.matrix,Rules=R,act.fun="sigmoidal",alpha=param[,2], beta=param[,3], method="lsoda")
