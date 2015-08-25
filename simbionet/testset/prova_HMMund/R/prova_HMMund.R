library(HMMb)

init()

MOD <- createMOD(m=2,auto=TRUE)
N = 5
system.time(ris <- HMMund(N=N,MODULES=MOD))
ris