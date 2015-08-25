library("HMM")

MOD <- createMOD(m=2,auto=TRUE)
system.time(HMMund(N=5, MODULES=MOD, sepgraph=FALSE))
cat("ok")