library("HMM")

#~ HMM<-function(N=50, Cf.cl=0.3, gamma=2.2, DEGREE=NULL, INdegree=c("free", "out"), MODULES, prior.p.subnet=NULL, max.con=12, feedback=TRUE, zero.nodes.with.indegree0=TRUE, sepgraph=TRUE, r.tol=0.1, a.tol=1, iter=1)

MOD <- createMOD(m=2,auto=TRUE)
p <- c(0.10417369, 0.25763653, 0.17851329, 0.19361568, 0.05220486, 0.10911545, 0.10474050)
system.time(HMMund(N=5, MODULES=MOD, prior.p.subnet=p))
cat("ok")