setwd("H:\\R\\r2c\\ok")

library(netsim)
dyn.load("probmod.dll")


probmod_C <-
    function(M, h, Sin, Sout, STin, STout, Freq.in, Freq.out, toll)
{
    stopifnot(is.matrix(M))
    stopifnot(is.vector(h))
    stopifnot(is.vector(Sin))
    stopifnot(is.vector(Sout))
    stopifnot(is.vector(STin))
    stopifnot(is.vector(STout))
    stopifnot(is.vector(Freq.in))
    stopifnot(is.vector(Freq.out))
    stopifnot(is.vector(toll))
    lista <- .Call("probmod",
                 as.matrix(M),
                 as.vector(h),
                 as.vector(Sin),
                 as.vector(Sout),
                 as.vector(STin),
                 as.vector(STout),
                 as.vector(Freq.in),
                 as.vector(Freq.out),
                 as.vector(toll)
					)
    return(lista)
}

dm = round(runif(1) * 10) + 1

m = matrix(0, dm, dm)
for (i in 1:dm) {
    for (j in 1:dm)
        m[i, j] = sample(c(0, 1), 1)
}
#~ m

h = rep(0, dm)
for (i in 1:dm)
    h[i] = sample(seq(1, dm), 1, replace=TRUE)
#~ h

Sin = rep(0, dm)
for (i in 1:dm)
    Sin[i] = runif(1)
#~ Sin

Sout = rep(0, dm)
for (i in 1:dm)
    Sout[i] = runif(1)
#~ Sout

STin = rep(0, dm)
for (i in 1:dm)
    STin[i] = runif(1)
#~ STin

STout = rep(0, dm)
for (i in 1:dm)
    STout[i] = runif(1)
#~ STout

Freq.in = rep(0, dm)
for (i in 1:dm)
    Freq.in[i] = runif(1)
#~ Freq.in

Freq.out = rep(0, dm)
for (i in 1:dm)
    Freq.out[i] = runif(1)
#~ Freq.out

toll = rep(0, dm)
for (i in 1:dm)
    toll[i] = runif(1)
#~ toll

l = probmod(m, h, Sin, Sout, STin, STout, Freq.in, Freq.out, toll)
l

l1 = probmod_C(m, h, Sin, Sout, STin, STout, Freq.in, Freq.out, toll)
l1

cat("Elementi diversi, pos. 1: ", which(l[[1]] != l1[[1]]), "\n")
cat("Elementi diversi, pos. 2: ", which(l[[2]] != l1[[2]]), "\n")
cat("Elementi diversi, pos. 3: ", which(l[[3]] != l1[[3]]), "\n")
cat("Elementi diversi, pos. 4: ", which(l[[4]] != l1[[4]]), "\n")