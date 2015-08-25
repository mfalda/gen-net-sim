setwd("H:\\R\\r2c\\ok")

library(netsim)
dyn.load("createRules.dll")

createRules_C <-
    function (matr, f_indx)
{
    stopifnot(is.matrix(matr))
    lista <- .Call("createRules",
                 as.matrix(matr),
                 as.integer(f_indx)
            )
    return(lista)
}

dm = sample(c(1,2,3,4), 1)

M = matrix(0, dm, dm)
for (i in 1:dm) {
    for (j in 1:dm)
        M[i, j] = runif(1)
}
M

l = createRules(M, NULL)
l

l1 = createRules_C(M, -1)
l1