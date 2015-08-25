library("HMM")

#~ HMMund<-function(N=50, Cf.cl=0.3, gamma=2.2, DEGREE=NULL, MODULES, prior.p.subnet=NULL, max.con=12, sepgraph=TRUE, r.tol=0.1, a.tol=1, iter=1)

MOD <- createMOD(m=2,auto=TRUE)
system.time(HMMund(N=5, MODULES=MOD))
cat("ok")