library(HMMb)

init()

MOD <- createMOD(m=2,auto=TRUE)
N = 5
system.time(ris <- HMM(N=N,MODULES=MOD, zero.nodes.with.indegree0=FALSE))
ris
