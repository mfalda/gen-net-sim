library(HMMb)

init(F, F)
MOD <- createMOD(m=2,auto=TRUE)
N <- 1000
L <- length(MOD)
p <- runif(L)
p <- p / sum(p)
deg <- c()
somma <- N
for (i in 1:N) {
    tmp <- runif(1, 0, (N + 1) / N)
    if (somma - tmp > 0) {
        somma <- somma - tmp
    }
    else tmp <- 0.0
    deg <- c(deg, tmp)
}
if (somma > 0) deg <- c(deg, somma) else deg <- c(deg, 0.0)
#cat("p = ", p)
#cat("deg = ", deg)
sum(deg)

fp <- file("stat_R.tsv", "w")
P <- 10
s <- 2 ^ seq(1, 9)

for (i in 1:P) {
    cat(file=fp, sprintf("\t%d", i))
}
cat(file=fp, "\tMedia\tDev. std.\n")

times <- c()
for (n in s) {
    cat(file=fp, sprintf("%d", n))
    for (i in 1:P) {
        t <- system.time(ris <- HMM(N=n,MODULES=MOD,Cf.cl=0, iter=1))
        cat(file=fp, sprintf("\t%3.3f", t[1]))
        times <- c(times, t[1])
    }
    cat(file=fp, sprintf("\t%3.3f\t%3.3f\n", mean(times), sd(times)))
}
close(fp)