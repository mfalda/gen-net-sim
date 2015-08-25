library(HMMb)

init()
leggi_prob()

MOD <- createMOD(m=2,auto=TRUE)
MOD
N = 5
system.time(ris <- HMM(N=N,MODULES=MOD, gamma=1.2))
ris
