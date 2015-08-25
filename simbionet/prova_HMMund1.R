library(HMM)

MOD <- createMOD(m=2,auto=TRUE)
N = 5
p = c(0.1055624, 0.1900999, 0.07080406, 0.1529689, 0.1150041, 0.04754141, 0.3180193)
dg = c(0.8803356, 0.5722373, 0.3792771, 0.4924523, 0.5471738, 2.128524)
cat("somma = ", sum(dg,na.rm=TRUE))
system.time(ris <- HMMund(N=N,MODULES=MOD,Cf.cl=0, DEGREE=dg, iter=1))
#~ system.time(ris <- HMMund(N=N,MODULES=MOD))
cat("ok")