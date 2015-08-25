library(HMM)

MOD <- createMOD(m=2,auto=TRUE)
N = 5
p <- runif(N + 2)
p <- p / sum(p)
#~ dg = c(0.8803356, 0.5722373, 0.3792771, 0.4924523, 0.5471738, 2.128524)
system.time(ris <- HMMund(N=N,MODULES=MOD,Cf.cl=0, prior.p.subnet=p, iter=1))
#~ system.time(ris <- HMMund(N=N,MODULES=MOD))
cat("ok")