library(HMMb)

init()

MOD <- createMOD(m=2,auto=TRUE)
N = 5
p <- c(0.10417369, 0.25763653, 0.17851329, 0.19361568, 0.05220486, 0.10911545, 0.10474050)
system.time(HMMund(N=5, MODULES=MOD, prior.p.subnet=p))
ris