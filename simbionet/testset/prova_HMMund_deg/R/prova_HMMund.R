library(HMMb)

init()

#~ HMM<-function(N=50, Cf.cl=0.3, gamma=2.2, DEGREE=NULL, INdegree=c("free", "out"), MODULES, prior.p.subnet=NULL, max.con=12, feedback=TRUE, zero.nodes.with.indegree0=TRUE, sepgraph=TRUE, r.tol=0.1, a.tol=1, iter=1)

MOD <- createMOD(m=2,auto=TRUE)
N <- 5
deg <- c(0.6759626, 0.8962194, 0.2949845, 0.1538275, 1.1431788, 1.8358272)
system.time(HMMund(N=5, MODULES=MOD, DEGREE=deg))
cat("ok")
