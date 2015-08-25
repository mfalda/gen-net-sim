setwd("H:\\R\\r2c\\ok")

library(netsim)
dyn.load("probmod.dll")

Score_C <-
    function(S, ST, Freq, n, toll)
{
    stopifnot(is.vector(S))
    stopifnot(is.vector(ST))
    stopifnot(is.vector(Freq))
    #~ stopifnot(is.integer(n))
    stopifnot(is.vector(toll))
    vett <- .Call("Score",
                 as.vector(S),
                 as.vector(ST),
                 as.vector(Freq),
                 as.integer(n),
                 as.vector(toll)
					)
    return(vett)
}

dm = round(runif(1) * 10) + 1
S = rep(0, dm)
for (i in 1:dm)
    S[i] = sample(seq(0, dm), 1)
S

ST = rep(0, dm)
for (i in 1:dm)
    ST[i] = runif(1)
ST

Freq = rep(0, dm)
for (i in 1:dm)
    Freq[i] = runif(1)
Freq

n = dm#round(runif(1) * 10) + 1
n

toll = rep(0, dm)
for (i in 1:dm)
    toll[i] = runif(1)
toll

ris = Score(S, ST, Freq, n, toll)
ris

ris1 = Score_C(S, ST, Freq, n, toll)
ris1

cat("Elementi diversi: ", which(ris != ris1), "\n")

