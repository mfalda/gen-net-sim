library(HMM)

MOD <- createMOD(m=2,auto=TRUE)
N = 5
deg <- c()
somma <- N
for (i in 1:N) {
    if (somma > 0) {
        tmp <- runif(1, 0, (N + 1) / N)
        somma <- somma - tmp
    }
    else tmp <- 0.0
    deg <- c(deg, tmp)
}
if (somma > 0) deg <- c(deg, somma) else deg <- c(deg, 0.0)
system.time(ris <- HMM(N=N,MODULES=MOD,Cf.cl=0, DEGREE=deg, iter=1))
#~ system.time(ris <- HMMund(N=N,MODULES=MOD))
cat("ok")