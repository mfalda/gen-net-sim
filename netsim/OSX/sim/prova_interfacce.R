setwd("/Users/mf/Desktop/linux/sim")

dyn.load("interfacce.so")

ret_vettore <-
    function (v, num)
{
    stopifnot(is.vector(v))
    stopifnot(is.numeric(v))
    vett <- .Call("ret_vettore",
                 as.vector(v),
                as.double(num)
                 )
    return(vett)
}

ret_matrice <-
    function (m)
{
    stopifnot(is.matrix(m))
    matr <- .Call("ret_matrice",
                 as.matrix(m)
                 )
    return(matr)
}

ret_lista <-
    function (l)
{
    stopifnot(is.list(l))
    lst <- .Call("ret_lista",
                 as.list(l)
                 )
    return(lst)
}

v = c(1, 2, 3)
v

ris <- ret_vettore(v, 10.5)
ris

m = matrix(c(-1, -2, 3, 10, -4, 1), 2, 3)
m

ris <- ret_matrice(m)
ris

v1 <- c(11.2, 23.45, 32.32)
l = list()
l[[1]] = v1
l[[2]] = m
l[[3]] = "stringa"
l

ris <- ret_lista(l)
ris
