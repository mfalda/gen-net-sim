setwd("H:\\R\\r2c\\ok")


library("netsim")
dyn.load("cluster_coeff.dll")

cluster.coeff_C <-
    function(W)
{
    stopifnot(is.matrix(W))
    lista <- .Call("cluster_coeff",
                 as.matrix(W)
					)
    return(lista)
}

dm = round(runif(1) * 10) + 1

w = matrix(0, dm, dm)
for (i in 1:dm) {
    for (j in 1:dm)
        w[i, j] = runif(1)
}
w

l1 = cluster.coeff_C(w)
l1

l = cluster.coeff(w)
l

