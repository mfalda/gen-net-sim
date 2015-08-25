library(HMMb)

init()

MOD <- createMOD(m=2,auto=TRUE)
MOD
N = 5
system.time(ris <- HMM(N=N,MODULES=MOD, Cf.cl=0.0))
ris
